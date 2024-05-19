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
    std::vector<double> file;
    int temp_a = a.getNumChannels();
    int temp_b = a.getNumSamplesPerChannel();
    for (int i = 0; i < temp_a; i++)
    {
        for (int j = 0; j < temp_b; j++)
        {
            file.push_back(a.samples[i][j]);
        }
    }
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
    static int nr_pliku=0;
    plot(X, Y);
    save("plot"+std::to_string(nr_pliku)+".png");
    nr_pliku++;
}




std::vector<std::complex<double>> DFT(std::vector<double> input )
{
    
    std::vector<std::complex<double>> wynik4;
    std::complex<double> sum;
    int N = input.size();
    int k = N;

    for (int K = 0; K < k; K++) {
        sum.real(0);
        sum.imag(0);
        for (int n = 0; n < N; n++)
        {
            std::complex<double> temp;
            temp.real(input[n] * cos((2.f * M_PI * K * static_cast<float> (n)) / N));
            temp.imag(input[n] * -sin((2.f * M_PI * K * static_cast<float> (n)) / N));
            sum += temp;
        }
        wynik4.push_back(sum);
    }
    return  wynik4;
}

std::vector<double> IDFT(std::vector<std::complex<double>> input) 
{
    std::complex<double> sum;
    int N = input.size();
    std::vector<double> inverse(N);

    for (int n = 0; n < N; n++) {

        sum.real(0);
        sum.imag(0);
        for (int k = 0; k < N; k++) {
            double wnetrze = 2.f * M_PI * k * n / N; //e^{ix}=cos{x}+isin{x}
            std::complex<double> euler(cos(wnetrze), sin(wnetrze)); 
            sum += input[k] * euler; // an*Wn^(...), Wn=e^i*(2pi/N)
        }
        inverse[n] = sum.real() / N;
    }

    return inverse;
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
    module.def("IDFT", &IDFT);
    module.def("generate_sine_wave", &generate_sine_wave);
}