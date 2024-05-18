#define _USE_MATH_DEFINES
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <stdio.h>
#include <matplot/matplot.h>
#include <cmath>
#include <string>
#include <numbers>
#include <AudioFile.h>

namespace py = pybind11;

void say_hello(int x) { //funkcja do testowania czy biblioteka dziala

    printf("Hello World!\n");
    std::cout << x << std::endl;

}

void writeSineWaveToAudioFile(float sampleRate = 44100.f , float frequencyInHz = 440.f)
{
    AudioFile<float> a;
    a.setNumChannels(2);
    a.setNumSamplesPerChannel(44100);


    // 3. Write the samples to the AudioFile sample buffer

    for (int i = 0; i < a.getNumSamplesPerChannel(); i++)
    {
        for (int channel = 0; channel < a.getNumChannels(); channel++)
        {
            a.samples[channel][i] = sin((static_cast<float> (i) / sampleRate) * frequencyInHz * 2.f * (float)M_PI);
        }
    }

    // 4. Save the AudioFile

    std::string filePath = "sine-wave.wav"; // change this to somewhere useful for you
    a.save("sine-wave.wav", AudioFileFormat::Wave);
}

void writeCosineWaveToAudioFile(float sampleRate = 44100.f, float frequencyInHz = 440.f)
{
    AudioFile<float> a;
    a.setNumChannels(2);
    a.setNumSamplesPerChannel(44100);


    // 3. Write the samples to the AudioFile sample buffer

    for (int i = 0; i < a.getNumSamplesPerChannel(); i++)
    {
        for (int channel = 0; channel < a.getNumChannels(); channel++)
        {
            a.samples[channel][i] = sin((static_cast<float> (i) / sampleRate) * frequencyInHz * 2.f * (float)M_PI + (float)M_PI/2.f);
        }
    }

    // 4. Save the AudioFile

    std::string filePath = "cosine-wave.wav"; // change this to somewhere useful for you
    a.save("cosine-wave.wav", AudioFileFormat::Wave);
}

void writeSquareWaveToAudioFile(float sampleRate = 44100.f, float frequencyInHz = 440.f)
{
    AudioFile<float> a;
    a.setNumChannels(2);
    a.setNumSamplesPerChannel(sampleRate);


    // 3. Write the samples to the AudioFile sample buffer

    for (int i = 0; i < a.getNumSamplesPerChannel(); i++)
    {
        for (int channel = 0; channel < a.getNumChannels(); channel++)
        {
            if (sin((static_cast<float> (i) / sampleRate) * frequencyInHz * 2.f * (float)M_PI) > 0) {
                a.samples[channel][i] = 1.f;
            }else
                if (sin((static_cast<float> (i) / sampleRate) * frequencyInHz * 2.f * (float)M_PI) <= 0) {
                    a.samples[channel][i] = -1.f;
                }
                
            
        }
    }

    // 4. Save the AudioFile

    std::string filePath = "square-wave.wav"; // change this to somewhere useful for you
    a.save("square-wave.wav", AudioFileFormat::Wave);
}


void loadAudioFileAndProcessSamples(std::string inputFilePath )
{
    
   

    //---------------------------------------------------------------
    // 2. Create an AudioFile object and load the audio file

    AudioFile<float> a;
    bool loadedOK = a.load(inputFilePath);

    /** If you hit this assert then the file path above
     probably doesn't refer to a valid audio file */
    assert(loadedOK);

    //a.samples[channel][i] = 
   
    

    using namespace matplot;

    plot(a.samples);

    show();

}

void generate_sine_wave(std::vector <double> wynik ,double frequency=1,int sample_rate=10000) {

    std::vector <double> wynik1;
        for(int i=0;i<sample_rate;i++)
            wynik1.push_back(sin((static_cast<float> (i) / sample_rate) * frequency * 2.f * (float)M_PI)); 
}

void generate_cosine_wave(std::vector <double> wynik, double frequency , int sample_rate ) {

    std::vector <double> wynik2;
    for (int i = 0; i < sample_rate; i++)
        wynik2.push_back(cos((static_cast<float> (i) / sample_rate) * frequency * 2.f * (float)M_PI));
}

void generate_square_wave(std::vector <double> wynik, double frequency , int sample_rate ) {

    std::vector <double> wynik3;
    for (int i = 0; i < sample_rate; i++)
        if (sin((static_cast<float> (i) / sample_rate) * frequency * 2.f * (float)M_PI) > 0) {
            wynik3.push_back(1);
        }
        else
            if (sin((static_cast<float> (i) / sample_rate) * frequency * 2.f * (float)M_PI) <= 0) {
                wynik3.push_back(-1) ;
            }
        
}

void plot_line(std::vector<double> X, std::vector<double> Y) {

    using namespace matplot;
    
    plot(X, Y);
    save("plot.png");
}




std::vector<std::complex<double>> DFT(std::vector<double> input )
{
    double a = 0;
    double b = 0;
    std::vector<std::complex<double>> wynik4;
    int N = input.size();
    float k = N;
    for (int n = 0; n < N; n++)
    {
        wynik4.push_back(0);
        wynik4[n].real(input[n] * cos((2.f * M_PI * k * static_cast<float> (n)) / N))  ;
        wynik4[n].imag( input[n] * -sin((2.f * M_PI * k * static_cast<float> (n)) / N)) ;
    }
    std::complex<double> temp(a, b);
    return  wynik4;
}



PYBIND11_MODULE(pybind11module, module) {

    module.doc() = "Pybind11Module";

    module.def("say_hello", &say_hello);
    module.def("sin", &writeSineWaveToAudioFile);
    module.def("cos", &writeCosineWaveToAudioFile);
    module.def("square", &writeSquareWaveToAudioFile);
    module.def("test", &loadAudioFileAndProcessSamples);
    module.def("plot_line", &plot_line);
    module.def("DFT", &DFT);
    module.def("generate_sine_wave", &generate_sine_wave);
}