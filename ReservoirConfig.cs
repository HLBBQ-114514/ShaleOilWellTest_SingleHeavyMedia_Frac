using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ShaleOilWellTest
{
    internal class ReservoirConfig(double h, double k, double Bo, double mu,
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
        public double p_D { get; set; }
        public double t_D { get; set; }
        public double h_t { get; set; }
        public double omega { get; set; }
        public double avgLambda { get; set; }
        public double avgphiCt { get; set; }
        public double eta { get; set;}
        public double zeta { get; set; }
        public double avgeta { get; set; }
        public double CsD { get; set; }
        public double reD { get; set; }
        public void factorlessness(ReservoirConfig config)
        {
            double Lref = rw;
             reD = re / Lref;
             h_t = h + config.h;
            //平均流度
            avgLambda = h / h_t * (config.k / config.mu + k / mu);
            //平均储容系数
            avgphiCt = 1 / h_t * (phi * ct * h + config.phi * config.ct * config.h);
            //平均导压系数
            eta = (avgLambda/ avgphiCt);

            zeta = k / mu * h / (avgLambda * h_t);
            omega = phi * ct * h / (avgphiCt * h_t);

            avgeta = zeta / omega;
            CsD = Cs / (6.2832 * avgphiCt * h_t * Lref * Lref);
        }
    }
}
