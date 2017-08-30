// psk_orange
// A homebrew psk31 codec
// by John howard, orthopteroid@gmail.com, August 2017
// MIT license

// for additional info see:
// https://en.wikipedia.org/wiki/PSK31
// http://aintel.bi.ehu.es/psk31.html

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <deque>
#include <vector>
#include <functional>
#include <cstring>

#include "PSKOrange.h"

using namespace std;

//////////////////////////
// int16 LE file writing

void write_s16LE(ostream &stream, const float &value)
{
    int16_t s16 = (uint16_t)((float)INT16_MAX * value);
    stream.write(reinterpret_cast<const char*>(&s16), sizeof(int16_t));
}

float read_s16LE(istream &stream)
{
    int16_t s16;
    stream.read(reinterpret_cast<char*>(&s16), sizeof(int16_t));
    return (float)s16 / (float)INT16_MAX;
}

//////////////////////////

int main()
{
    srand48(int(time(NULL)));

    const int element_cycles = 50;
    const int carrier_cycles = 3 * element_cycles;
    const int freq = 500;
    const int samplerate = 22050;

    // basic test
    if(false)
    {
        // synthesis
        {
            ofstream outFile;
            outFile.open("out.bin", ios::binary);

            bool q[] = {1,0,1,0,1,0};

            PSKOrange<element_cycles,carrier_cycles,freq,samplerate>::ElementEncoder encoder;
            encoder.encode(
                [&outFile](float v) { write_s16LE( outFile, v ); },
                [&q](void) -> bool* { static int i = 0; return i < 6 ? &q[i++] : nullptr; }
            );

            outFile.close();

            //system("/usr/bin/aplay -f cd -r 22050 out.bin"); // enable to play audio!
        }

        // read and decode
        {
            ifstream inFile;
            inFile.open("out.bin", ios::binary);

            std::deque<bool> bitseq;
            std::ostream &console = cout; // lower cout into local context for use in lambda

            PSKOrange<element_cycles,carrier_cycles,freq,samplerate>::ElementDecoder decoder;
            decoder.decode(
                [&inFile](void) -> float*
                {
                    static float sample;
                    return inFile.eof() ? nullptr : &(sample = read_s16LE(inFile));
                },
                [&inFile,&bitseq](bool bit) -> void { bitseq.push_back(bit); },
                [&console](int bit) -> void { /* static int x = 0; console << "error in bit " << bit << " " << x++ << std::endl; */ } // todo: fix
            );

            // drop last bit as it is the transition from last element to the ending carrier
            bitseq.pop_back();

            for(auto iter = bitseq.begin(); iter != bitseq.end(); iter++) { cout << *iter << " "; }
            cout << std::endl;

            inFile.close();
        }
    }

    // full message string
    if(true)
    {
        // synthesis
        {
            ofstream outFile;
            outFile.open("out.bin", ios::binary);

            const char* message = "Hello World!";

            std::ostream &console = cout; // lower cout into local context for use in lambda

            PSKOrange<element_cycles,carrier_cycles,freq,samplerate>::ElementEncoder encoder;
            encoder.encode(
                [&outFile](float v) { write_s16LE( outFile, v ); },
                [&message,&console](void) -> bool*
                {
                    static int i = 0;
                    static bool bit = false;

                    if(message[i/8] == 0) { return nullptr; }
                    bit = (bool)(message[i/8] & (1 << (i%8)));

                    console << bit << std::flush;
                    if(++i%8 == 0) { console << ' ' << std::flush; }

                    return &bit;
                }
            );

            outFile.close();

            //system("/usr/bin/aplay -f cd -r 22050 out.bin"); // enable to play audio!
        }
        cout << std::endl;

        // read and decode
        {
            ifstream inFile;
            inFile.open("out.bin", ios::binary);

            std::ostream &console = cout; // lower cout into local context for use in lambda

            PSKOrange<element_cycles,carrier_cycles,freq,samplerate>::ElementDecoder decoder;
            decoder.decode(
                    [&inFile](void) -> float*
                    {
                        static float sample; // simple hack that allows us to signal end-of-stream
                        return inFile.eof() ? nullptr : &(sample = read_s16LE(inFile));
                    },
                    [&console](bool bit) -> void
                    {
                        static int i = 0;

                        console << bit << std::flush;
                        if(++i%8 == 0) { console << ' ' << std::flush; }
                    },
                    [&console](int bit) -> void { /*console << "(E" << bit << ')';*/ } // todo: fix
            );

            inFile.close();
        }
        cout << std::endl;
    }

    // read and generate details
    if(true)
    {
        ifstream inFile;
        inFile.open("out.bin", ios::binary);

        ofstream detailsFile;
        detailsFile.open("details.bin", ios::binary);

        PSKOrange<element_cycles,carrier_cycles,freq,samplerate>::ElementDecoder decoder;
        while(!inFile.eof())
        {
            float sample = read_s16LE(inFile);
            decoder.process(sample);

            write_s16LE(detailsFile, sample);
            write_s16LE(detailsFile, decoder.inversionDetector.phase * M_1_PI);

            write_s16LE(detailsFile, decoder.channelA.phase * M_1_PI);
            write_s16LE(detailsFile, decoder.channelA.iIntegrator);
            write_s16LE(detailsFile, decoder.channelA.qIntegrator);
            write_s16LE(detailsFile, decoder.channelA.squelch);

            write_s16LE(detailsFile, decoder.channelB.phase * M_1_PI);
            write_s16LE(detailsFile, decoder.channelB.iIntegrator);
            write_s16LE(detailsFile, decoder.channelB.qIntegrator);
            write_s16LE(detailsFile, decoder.channelB.squelch);
        }

        inFile.close();
        detailsFile.close();
    }

    return 0;
}
