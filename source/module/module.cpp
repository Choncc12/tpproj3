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

std::vector<double> loadAudioFileAndProcessSamples(std::string inputFilePath, int length)
{
    
    AudioFile<float> a;
    bool loadedOK = a.load(inputFilePath);

    assert(loadedOK);
    std::vector<double> file;
    int temp_a = a.getNumChannels();
    int temp_b = a.getNumSamplesPerChannel();
    for (int i = 0; i < temp_a; i++)
    {
        for (int j = 0; j < temp_b; j++)
        {
            file.push_back(a.samples[i][j]);
            if (file.size() >= length) {
                return file;
            }
        }
    }
    return file;
}

std::vector <double> generate_sine_wave(double frequency = 1, int sample_rate = 10000) {

    std::vector <double> wynik1;
    for (int i = 0; i < sample_rate; i++) {
        wynik1.push_back(sin((static_cast<float> (i) / sample_rate) * frequency * 2.f * (float)M_PI));
    }
    return wynik1;
}

std::vector <double> generate_cosine_wave(double frequency = 1, int sample_rate = 10000) {

    std::vector <double> wynik2;
    for (int i = 0; i < sample_rate; i++) {
        wynik2.push_back(cos((static_cast<float> (i) / sample_rate) * frequency * 2.f * (float)M_PI));
    }
    return wynik2;
}

std::vector <double> generate_square_wave(double frequency = 1, int sample_rate = 10000) {

    std::vector <double> wynik3;
    for (int i = 0; i < sample_rate; i++) {
        if (sin((static_cast<float> (i) / sample_rate) * frequency * 2.f * (float)M_PI) >= 0) {
            wynik3.push_back(1);
        }
        else {
            wynik3.push_back(-1);
        }
    }
    return wynik3;
}

std::vector <double> generate_sawtooth_wave(double frequency = 1, int sample_rate = 10000) {

    std::vector <double> wynik4;
    for (int i = 0; i < sample_rate; i++) {
        wynik4.push_back(fmod(static_cast<float>(i), frequency));
    }
    return wynik4;
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

std::vector <double> low_pass_filter(std::vector <double>input, double frequency) {

    std::vector<std::complex<double>> w_dft;
    w_dft = DFT(input);
    int b = ceil(frequency);
    for (int i = b; i < w_dft.size(); i++) {
        w_dft[i] = 0;
    }
    std::vector<double> w_idft;
    w_idft = IDFT(w_dft);
    return w_idft;
}

PYBIND11_MODULE(pybind11module, module) {

    module.doc() = "Pybind11Module";

    module.def("audio", &loadAudioFileAndProcessSamples);
    module.def("plot_line", &plot_line);
    module.def("DFT", &DFT);
    module.def("IDFT", &IDFT);
    module.def("generate_sine_wave", &generate_sine_wave);
    module.def("generate_cosine_wave", &generate_cosine_wave);
    module.def("generate_square_wave", &generate_square_wave);
    module.def("generate_sawtooth_wave", &generate_sawtooth_wave);
    module.def("low_pass_filter", &low_pass_filter);
}