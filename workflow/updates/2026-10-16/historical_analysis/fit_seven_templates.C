// Seven-channel adaptation of Bill's fun_sim_total histogram-template model.
// Source: https://github.com/billlee77/omega_analysis/blob/master/fitting_functions.C
// No interpolation, peak shifts, width changes, or additional background.
#include <TH1.h>
TH1* seven_templates[7] = {nullptr};
void set_seven_template(int j, TH1* h) { seven_templates[j] = h; }
double seven_template_total(double* x, double* p) {
    double result = 0;
    for (int j = 0; j < 7; ++j) {
        int bin = seven_templates[j]->GetXaxis()->FindFixBin(x[0]);
        result += p[j] * seven_templates[j]->GetBinContent(bin);
    }
    return result;
}
