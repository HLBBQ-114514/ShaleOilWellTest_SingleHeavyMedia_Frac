using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;

namespace ShaleOilWellTest
{
    internal class ReservoirConfigMedia(double h, double k, double Bo, double mu,
        double ct, double Cs, double rw, double re, double phi, double s, double sf,double xf)
    {
        public static double Pi;
        public double h = h;
        public double k = k;
        public double Bo = Bo;
        public double mu = mu;
        public double ct = ct;
        public double Cs = Cs;
        public double rw = rw;
        public double re = re;
        public double phi = phi;
        public double s = s;
        public double sf = sf;
        public double xf = xf;

        //无因次化
        public double h_t { get; set; }
        public double omega { get; set; }
        public double avgLambda { get; set; }
        public double avgphiCt { get; set; }
        public double eta { get; set;}
        public double zeta { get; set; }
        public double avgeta { get; set; }
        public double CsD { get; set; }
        public double reD { get; set; }
        public void factorlessness(ReservoirConfigMedia config,double omega,double zeta)
        {
             reD = re / xf;
             h_t = h + config.h;
            //平均流度
            avgLambda = h / h_t * (config.k / config.mu + k / mu);
            //平均储容系数
            avgphiCt = 1 / h_t * (phi * ct * h + config.phi * config.ct * config.h);
            //平均导压系数
            eta = (avgLambda/ avgphiCt);

            this.zeta = zeta;// k / mu * h / (avgLambda * h_t);
            this.omega = omega;//phi * ct * h / (avgphiCt * h_t);

            avgeta = zeta / omega;
            CsD = Cs / (6.2832 * avgphiCt * h_t * xf * xf);
        }
    }

    internal class ReservoirConfigDoubleMedia(double h, double k, double Bo, double mu,
        double ctm,double ctf, double Cs, double rw, double re, double phiM, double phiF, double s, double sf, double xf)
    {
        public static double Pi;
        public double h = h;
        public double kf = k;
        public double Bo = Bo;
        public double mu = mu;
        public double ctm = ctm;
        public double ctf = ctf;
        public double Cs = Cs;
        public double rw = rw;
        public double re = re;
        public double phiM = phiM;
        public double phiF = phiF;
        public double s = s;
        public double sf = sf;
        public double xf = xf;
        public double theta { get; set; }
        public double xfD;
        public double h_t { get; set; }
        public double reD { get { return re / rw; } set { } }
        public double avgLambda { get; set; }
        public double avgphiCt { get; set; }
        public double eta { get; set; }
        public double omega { get; set; }
        public double zeta { get; set; }
        public double totalPhiCt { get; set; }
        public double omegaM=0.95;    
        public double omegaF=0.05;     
        public double lambdaF { get; set; }
        public void factorlessness(ReservoirConfigDoubleMedia[] configs)
        {
            double phiCtM = 0.0;
            var baselayer = configs[0];
            reD = re / baselayer.rw;
            xfD = xf / baselayer.xf;
            for (int i = 0; i < configs.Length; i++)
            {
                h_t = configs[i].h;
                //平均流度
                avgLambda += (configs[i].kf / configs[i].mu) * configs[i].h;
                //平均储容系数
                avgphiCt += (configs[i].phiF * configs[i].ctf * configs[i].h);
                phiCtM += (configs[i].phiM * configs[i].ctm);
            }
            avgLambda *= 1 / h_t;
            avgphiCt *= 1 / h_t;
            //平均导压系数
            eta = (avgLambda / avgphiCt);
/*            this.omegaM = this.phiM * this.ctm / (avgphiCt * h_t + phiCtM);
            this.omegaF = 1 - this.omegaM;*/
            zeta = k / mu * h / (avgLambda * h_t);
            omega = phiF * ctf * h / (avgphiCt * h_t);
            Debug.WriteLine($"zeta:{ zeta}, omega:{ omega}");
        }
    }
}
