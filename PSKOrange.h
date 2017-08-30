// psk_orange: A homebrew psk31-style codec
// orthopteroid@gmail.com, August 2017, MIT license

// for additional info see:
// https://en.wikipedia.org/wiki/PSK31
// http://aintel.bi.ehu.es/psk31.html

#ifndef PSK_ORANGE_PSKORANGE_H
#define PSK_ORANGE_PSKORANGE_H

#include <stdint.h>
#include <cmath>
#include <deque>
#include <vector>
#include <functional>

using namespace std;

//////////////

template<int element_cycles, int carrier_cycles, int freq, int samplerate>
struct PSKOrange
{
    // pskorange uses a delta-bitcode to send bits. deltas are signalled
    // as phase changes from the previous bit. a starting bit of 0 is assumed,
    // which is represented as a fixed-duration carrier establishing a
    // reference phase angle.
    static constexpr float signal_phase = (float)M_PI;

    struct ElementEncoder
    {
        static constexpr float amplification = .95; // to be safe from saturation

        struct Fader
        {
            float theta = 0.f, delta = 0.f;

            Fader() = default;

            void fadein(int samples) { theta = 0.f; delta = (float)M_PI_2 / (float)samples; }
            void fadeout(int samples) { theta = (float)M_PI_2; delta = (float)M_PI_2 / (float)samples; }

            float operator()() {
                float z = sin(theta);
                theta += delta;
                return z;
            }
        } fader;

        struct Osc
        {
            float theta = 0.f;

            static constexpr float delta = (2.f * (float)M_PI * (float)freq) / (float)samplerate;

            Osc() = default;

            void reset(float phase) { theta=phase; }

            float operator()() {
                float z = sin(theta);
                theta += delta;
                return z;
            }
        } osc;

        //////////

        ElementEncoder() = default;

        // send a psk element of the specified phase that lasts for the specified number of cycles
        template<int cycles_total>
        void transmit(
            const std::function<void(float)> &senderfn,
            bool signal
        )
        {
            static constexpr int cycles_fade = 10;
            static constexpr int cycles_body = cycles_total - 20;
            static constexpr int samples_fade = (int)((float)cycles_fade * (float)samplerate / (float)freq);
            static constexpr int samples_body = (int)((float)cycles_body * (float)samplerate / (float)freq);

            osc.reset(signal ? PSKOrange::signal_phase : 0);

            fader.fadein(samples_fade);
            for(int i=0; i<samples_fade; i++) { senderfn( amplification * fader() * osc() ); }

            for(int i=0; i<samples_body; i++) { senderfn( amplification * osc() ); }

            fader.fadeout(samples_fade);
            for(int i=0; i<samples_fade; i++) { senderfn( amplification * fader() * osc() ); }
        }

        // converts the bitsequence into a series of bit-inversions which are transmitted as a series of phase-inversions
        void encode(
            const std::function<void(float)> &senderfn,
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
        static constexpr int ringSize = (int)(4.f * (float)samplerate / (float)freq);

        // psk_orange de/encodes bit transitions. This is the value of the last bit decoded.
        bool last_decoded_bit = false;

        // a ring buffer sized to fit a single wavelength of the element frequency
        struct RingBuffer
        {
            std::vector<float> ring;
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
        struct PhaseDetector
        {
            float iIntegrator = 0.f, qIntegrator = 0.f; // state variables, for diagnostics
            float phase = 0.f;
            float squelch = 0.f;

            // compute indices for where the quadrants are in the waveform
            static constexpr int q1End = (int)((float)samplerate / (float)freq / 4.f);
            static constexpr int q2End = (int)((float)samplerate / (float)freq / 2.f);
            static constexpr int q3End = (int)(3.f * (float)samplerate / (float)freq / 4.f);
            static constexpr int q4End = (int)((float)samplerate / (float)freq);

            // each channel scans half the ring
            static constexpr int channelSize = ringSize / 2;

            RingBuffer &sampler;

            explicit PhaseDetector(RingBuffer &s) : sampler(s) {}

            void process()
            {
                iIntegrator = qIntegrator = 0.f;

                // integrate and normalize
                for(int i=0; i<channelSize; i++) {
                    float sample = sampler.ring[ sampler.readIndex ];
                    if(++sampler.readIndex == ringSize) { sampler.readIndex = 0; }

                    iIntegrator += (sampler.theta < q2End) ? sample : -sample;
                    qIntegrator += (sampler.theta >= q1End) && (sampler.theta < q3End) ? sample : -sample;
                    if(++sampler.theta == q4End) { sampler.theta = 0; }
                }
                iIntegrator /= 1.f * (float)channelSize;
                qIntegrator /= 1.f * (float)channelSize;

                phase = quick_atan2( qIntegrator, iIntegrator );
                squelch = sqrt( iIntegrator * iIntegrator + qIntegrator * qIntegrator );
            }

            bool betterThan(PhaseDetector& other) { return squelch > 2. * other.squelch; }

            // Approximates atan2(y, x) with maximum error of 0.1620 degrees
            // https://stackoverflow.com/a/14100975/968363 - with dbz and nan fix and renormalization to (-pi,+pi)
            float quick_atan2( float y, float x )
            {
                static const uint32_t sign_mask = 0x80000000;
                static const float b = 0.596227f;

                // Extract the sign bits
                uint32_t ux_s  = sign_mask & (uint32_t &)x;
                uint32_t uy_s  = sign_mask & (uint32_t &)y;

                // Determine the quadrant offset
                auto q = (float)( ( ~ux_s & uy_s ) >> 29 | ux_s >> 30 );

                // Calculate the arctangent in the first quadrant
                float bxy_a = fabsf( b * x * y );
                float num = bxy_a + y * y;
                float atan_1q =  num / (x * x + bxy_a + num + .0001f); // .0001 to fix dbz and nan issue

                // Translate it to the proper quadrant
                uint32_t uatan_2q = (ux_s ^ uy_s) | (uint32_t &)atan_1q;
                return (q + (float &)uatan_2q) * (float)M_PI_2 - (float)M_PI; // [0,4) to [0,2pi) to (-pi,+pi)
            }

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

            static constexpr float noise_floor = signal_phase * .01f; // set a noise tolerance
            static constexpr float value_threshold = signal_phase * .4f; // set a threshold for 0/1-ness.

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

                    // because the signalling phase angle is pi and atan2 is not exact the calculated
                    // phase angle may oscillate between quadrant 2 and 3 - meaning an oscillation of +/- pi.
                    // so, use the inner fabs to stay in the upper quadrants.
                    float phase_delta_abs = fabsf(fabsf(phase_e) - fabsf(phase_c));

                    // smooth the maximums of the phase differences
                    phase = max(phase, phase_delta_abs) * .4f + phase * .6f;
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
                carrierDetector(channelA, channelB, carrier_cycles * (float)samplerate / (float)freq),
                inversionDetector(channelA, channelB, carrierDetector)
        { }

        void process(float value)
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
            const std::function<float*(void)> &readfn,
            const std::function<void(bool)> &bitfn,
            const std::function<void(int)> &noisefn)
        {
            int bitnum = 0;
            float* pSample = nullptr;

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
