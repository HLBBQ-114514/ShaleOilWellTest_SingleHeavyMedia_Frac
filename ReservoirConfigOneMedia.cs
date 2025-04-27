using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;

namespace ShaleOilWellTest
{
    internal class ReservoirConfigOneMedia(double h, double k, double Bo, double mu,
        double ct, double Cs, double rw, double re, double phi, double s, double sf,double xf)
    {
        public double Pi;
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
        public void factorlessness(ReservoirConfigOneMedia config,double omega,double zeta)
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
        public double Pi;
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
        public double h_t { get; set; }
        public double reD { get { return re / xf; } set { } }
        public double avgLambda { get; set; }
        public double avgphiCt { get; set; }
        public double eta { get; set; }
        public double omega { get; set; }
        public double zeta { get; set; }
        public double totalPhiCt { get; set; }
        public double omegaM => omega * 0.95;    
        public double omegaF => omega * 0.05;     
        public double lambdaF { get; set; }
        public void factorlessness(ReservoirConfigDoubleMedia[] configs)
        {
            reD = re / xf;
            for (int i = 0; i < configs.Length; i++)
            {
                h_t = configs[i].h;
                //平均流度
                avgLambda += (configs[i].kf / configs[i].mu) * configs[i].h;
                //平均储容系数
                avgphiCt += (configs[i].phiF * configs[i].ctf * configs[i].h);

            }
            avgLambda *= 1 / h_t;
            avgphiCt *= 1 / h_t;
            //平均导压系数
            eta = (avgLambda / avgphiCt);

            /*zeta = k / mu * h / (avgLambda * h_t);
            omega = phi * ct * h / (avgphiCt * h_t);*/

        }
    }
}
