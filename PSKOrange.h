// psk_orange: A homebrew psk31-style codec
// orthopteroid@gmail.com, August 2017, MIT license

// for additional info see:
// https://en.wikipedia.org/wiki/PSK31
// http://aintel.bi.ehu.es/psk31.html

#ifndef PSK_ORANGE_PSKORANGE_H
#define PSK_ORANGE_PSKORANGE_H

#include <stdint.h>
#include <deque>
#include <vector>
#include <functional>

using namespace std;

//////////////

template<int element_cycles, int carrier_cycles, int freq, int samplerate, class M>
struct PSKOrange
{
    typedef typename M::Valuetype Valuetype;

    // pskorange uses a delta-bitcode to send bits. deltas are signalled
    // as phase changes from the previous bit. a starting bit of 0 is assumed,
    // which is represented as a fixed-duration carrier establishing a
    // reference phase angle.
    static constexpr float signal_phase = M::PI;

    struct ElementEncoder
    {
        struct Fader
        {
            float theta = 0.f, delta = 0.f;

            Fader() = default;

            void fadein(int samples) { theta = 0.f; delta = M::PI_2 / (float)samples; }
            void fadeout(int samples) { theta = M::PI_2; delta = M::PI_2 / (float)samples; }

            float operator()() {
                float z = M::sin(theta);
                theta += delta;
                return z;
            }
        } fader;

        struct Osc
        {
            float theta = 0.f;

            static constexpr float delta = (2.f * M::PI * (float)freq) / (float)samplerate;

            Osc() = default;

            void reset(float phase) { theta=phase; }

            float operator()() {
                float z = M::sin(theta);
                theta += delta;
                return z;
            }
        } osc;

        //////////

        ElementEncoder() = default;

        // send a psk element of the specified phase that lasts for the specified number of cycles
        template<int cycles_total>
        void transmit(
            const std::function<void(Valuetype)> &senderfn,
            bool signal
        )
        {
            static constexpr int cycles_fade = 10;
            static constexpr int cycles_body = cycles_total - 20;
            static constexpr int samples_fade = (cycles_fade * samplerate) / freq;
            static constexpr int samples_body = (cycles_body * samplerate) / freq;

            osc.reset(signal ? PSKOrange::signal_phase : 0);

            fader.fadein(samples_fade);
            for(int i=0; i<samples_fade; i++) { senderfn( fader() * osc() ); }

            for(int i=0; i<samples_body; i++) { senderfn( osc() ); }

            fader.fadeout(samples_fade);
            for(int i=0; i<samples_fade; i++) { senderfn( fader() * osc() ); }
        }

        // converts the bitsequence into a series of bit-inversions which are transmitted as a series of phase-inversions
        void encode(
            const std::function<void(Valuetype)> &senderfn,
            const std::function<bool*(void)> &bitfn
        )
        {
            transmit<carrier_cycles>(senderfn, false); // always starts with carrier phase

            bool invert = false;
            bool prevBit = false; // ie carrier phase angle
            bool* pBit = 0;
            while((pBit = bitfn()) != nullptr)
            {
                if(*pBit != prevBit) { invert = !invert; } // detect bit transition
                transmit<element_cycles>(senderfn, invert);
                prevBit = *pBit;
            }

            transmit<carrier_cycles>(senderfn, false); // always ends with carrier phase
        }
    };

    ////////////////

    struct ElementDecoder
    {
        // ring buffer size should be multiple waveform lengths to prevent small values
        // of (i,q) being calculated during element attenuation.
        static constexpr int ringSize = ((4 * samplerate) / freq);

        // psk_orange de/encodes bit transitions. This is the value of the last bit decoded.
        bool last_decoded_bit = false;

        // a ring buffer sized to fit a single wavelength of the element frequency
        struct RingBuffer
        {
            std::vector<Valuetype> ring;
            int writeIndex = 0; // destination index for new data
            int readIndex = 0; // index for reading data from buffer. each channel reads half the buffer. in sequence.
            int theta = 0; // phase angle of the current waveform sample, as a sample index.

            explicit RingBuffer() : ring(ringSize) {}

            void reset() { writeIndex = readIndex = theta = 0; }

            void process(float value)
            {
                ring[writeIndex] = value;
                if(++writeIndex == ringSize) { writeIndex = 0; }
                readIndex = writeIndex; // ie, the oldest written value
            }
        };

        // determines the phase angle of the waveform in a ringbuffer
        // ringbuffer is sized to fit an integer number of waveforms at current sample rate.
        // Comment from: https://stackoverflow.com/a/27546385/968363
        // Eqn N.20 of: http://kom.aau.dk/group/05gr506/report/node29.html
        // Other work: https://www.embedded.com/design/configurable-systems/4212086/DSP-Tricks--Frequency-demodulation-algorithms-
        struct PhaseDetector
        {
            float iIntegrator = 0.f, qIntegrator = 0.f; // state variables, for diagnostics
            float phase = 0.f;
            float squelch = 0.f;

            // compute indices for where the quadrants are in the waveform
            static constexpr int q4End = samplerate / freq;
            static constexpr int q1End = q4End / 4;
            static constexpr int q2End = 2 * q1End;
            static constexpr int q3End = 3 * q1End;

            // each channel scans half the ring
            static constexpr int channelSize = ringSize / 2;

            RingBuffer &sampler;

            explicit PhaseDetector(RingBuffer &s) : sampler(s) {}

            void process()
            {
                iIntegrator = qIntegrator = 0.f;

                // integrate and normalize
                for(int i=0; i<channelSize; i++) {
                    Valuetype sample = sampler.ring[ sampler.readIndex ];
                    if(++sampler.readIndex == ringSize) { sampler.readIndex = 0; }

                    iIntegrator += (sampler.theta < q2End) ? sample : -sample;
                    qIntegrator += (sampler.theta >= q1End) && (sampler.theta < q3End) ? sample : -sample;
                    if(++sampler.theta == q4End) { sampler.theta = 0; }
                }
                iIntegrator /= channelSize;
                qIntegrator /= channelSize;

                phase = M::atan2( qIntegrator, iIntegrator );
                squelch = M::sqrt( iIntegrator * iIntegrator + qIntegrator * qIntegrator );
            }

            bool betterThan(PhaseDetector& other) { return squelch > 2 * other.squelch; }

        };

        // identifies a carrier by it's length, then either locks the carrier to that channel (in the
        // case of the leading carrier of a bit-sequence) or releases it (in the case of the trailing
        // carrier of a bit-sequence)
        struct CarrierDetector {
            int channel_monitor = 0; // 0 == none, 1 == A, 2 == B, -1 == A falling edge, -2 == B falling edge
            int channel_lock = 0; // 0 == none, 1 == A, 2 == B
            int sample_count = 0;
            int samples_in_carrier = 0;

            PhaseDetector &channelA, &channelB;

            CarrierDetector(PhaseDetector &ca, PhaseDetector &cb, int s) :
                    channelA(ca), channelB(cb),
                    samples_in_carrier(s)
            {}

            void reset() { channel_monitor = channel_lock = sample_count = 0; }

            bool detected() { return channel_lock != 0; }

            void process()
            {
                sample_count++;

                // trailing-edge detection of stronger channel
                channel_monitor = channelA.betterThan(channelB) ? 1 : (channelB.betterThan(channelA) ? 2 : -channel_monitor);
                if(channel_monitor >= 0) { return; } // no edge yet

                // lock the carrier to this channel if there has been no other channel transitions
                // (ie an uninterrupted carrier of the correct length)
                if(sample_count >= samples_in_carrier)
                {
                    // clear or lock channel, depending if this is a leading or a trailing carrier
                    channel_lock = channel_lock > 0 ? 0 : -channel_monitor;
                }
                if(channel_lock > 0) { sample_count = 0; } // ignore squelches during carrier lock

                // clear the edge-detected state
                channel_monitor = 0;
            }
        };

        // detects phase inversions between a carrier channel and the data channel
        struct InversionDetector {
            float phase = 0.f;
            int channel = 0; // 0 == none, 1 == A, 2 == B, -1 == A falling edge, -2 == B falling edge
            int latch = -2; // -2 == undefined, -1 == NOISE, 0 == FALSE, 1 == TRUE

            static constexpr Valuetype noise_floor = signal_phase / 100; // set a noise tolerance
            static constexpr Valuetype value_threshold = (signal_phase / 5) * 2; // set a threshold for 0/1-ness.

            PhaseDetector &channelA, &channelB;
            CarrierDetector &carrierDetector;

            InversionDetector(PhaseDetector &ca, PhaseDetector &cb, CarrierDetector &cd) :
                    carrierDetector(cd),
                    channelA(ca), channelB(cb)
            {}

            void reset() { phase = 0.f; channel = 0; latch = -2; }

            void process()
            {
                if(!carrierDetector.detected()) { return; } // wait for carrier

                // discriminate stronger channel and it's trailing-edge
                channel = channelA.betterThan(channelB) ? 1 : (channelB.betterThan(channelA) ? 2 : -channel);

                if(channel > 0) // a channel is active
                {
                    float phase_c = carrierDetector.channel_lock == 1 ? channelA.phase : channelB.phase;
                    float phase_e = channel == 1 ? channelA.phase : channelB.phase;

                    // because the signalling phase angle is pi and atan2 is not exact, the calculated
                    // phase angle may oscillate between quadrant 2 and 3 - meaning an oscillation of +/- pi.
                    // so, we make liberal use of fabs to stay in the upper quadrants.
                    // Similar to Angle Unwrapping? See Fig N.4 of: http://kom.aau.dk/group/05gr506/report/node29.html
                    float phase_delta_abs = M::abs(M::abs(phase_e) - M::abs(phase_c));

                    // smooth the maximums of the phase differences
                    phase = M::max(phase, phase_delta_abs) * .4f + phase * .6f;
                }
                else if(channel < 0) // switching channels, determine bit value and latch it
                {
                    // if phase is above noise floor, determine signal value from phase angle.
                    // otherwise indicate a noise event.
                    latch = phase > noise_floor ? (phase > value_threshold ? +1 : 0) : -1;
                    channel = 0; // indicate no channel is active
                    phase = 0.f; // reset active phase
                }
            }

            // sample, clear and return latch value
            int read()
            { int sv = latch; latch = -2; return sv; }
        };

        RingBuffer ringBuffer;
        PhaseDetector channelA, channelB;
        CarrierDetector carrierDetector;
        InversionDetector inversionDetector;

        //////////

        ElementDecoder() :
                channelA(ringBuffer),
                channelB(ringBuffer),
                carrierDetector(channelA, channelB, (carrier_cycles * samplerate) / freq),
                inversionDetector(channelA, channelB, carrierDetector)
        { }

        void process(Valuetype value)
        {
            ringBuffer.process(value);
            channelB.process();
            channelA.process();
            carrierDetector.process();
            inversionDetector.process();
        }

        void reset()
        {
            ringBuffer.reset();
            inversionDetector.reset();
            carrierDetector.reset();
            last_decoded_bit = false; // corresponding to the state where the carrier is not yet detected
        }

        // return a decoded bit (0 or 1) or a value indicating no-value (-2) or noisy-signal (-1)
        int decode()
        {
            if(!carrierDetector.channel_lock) { return -2; } // invalid value when no carrier

            int invert = inversionDetector.read();
            if(invert < 0) { return invert; } // indicate invalid or noise

            if(invert == 1) last_decoded_bit = !last_decoded_bit; // decode delta
            return last_decoded_bit ? 1 : 0;
        }

        // read from a source until nullptr and send back any detected bits or noise events.
        void decode(
            const std::function<Valuetype*(void)> &readfn,
            const std::function<void(bool)> &bitfn,
            const std::function<void(int)> &noisefn)
        {
            int bitnum = 0;
            Valuetype* pSample = nullptr;

            while((pSample = readfn()) != nullptr)
            {
                process(*pSample);
                int value = decode();
                if(value >= 0) { bitfn((bool)value); bitnum++; }
                else if(value == -1) { noisefn(bitnum); }
            }
        }
    };

};

#endif //PSK_ORANGE_PSKORANGE_H
