using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ShaleOilWellTest
{
    internal class ReservoirConfig
    {
        public double Pi;
        public double h;
        public double k;
        public double Bo;
        public double mu;
        public double ct;
        public double Cs;
        public double rw;
        public double re;
        public double phi;
        public double s;
        public double sf;

        public ReservoirConfig(double h, double k, double Bo, double mu, 
            double ct, double Cs, double rw, double re, double phi, double s, double sf) 
        {
            this.h = h;
            this.k = k;
            this.Bo = Bo;
            this.mu = mu;
            this.ct = ct;
            this.Cs = Cs;
            this.rw = rw;
            this.phi = phi;
            this.s = s;
            this.re = re;
            this.sf = sf;
        }

        //无因次化
        public double p_D { get; set; }
        public double t_D { get; set; }
        public double h_t { get; set; }
        public double omega { get; set; }
        public double avgLambda { get; set; }
        public double avgphiCt { get; set; }
        public double eta { get; set;}
        public double zeta { get; set; }
        public void factorlessness(ReservoirConfig config)
        {   

            h_t = h + config.h;
            //平均流度
            avgLambda = h / h_t * (config.k/config.mu+k/mu);
            //平均储容系数
            avgphiCt = h/h_t*(phi*ct+config.phi*config.ct);
            //平均导压系数
            eta = (avgLambda/ avgphiCt);

            zeta = k / mu * h / (avgLambda * h_t);
            omega = phi * ct * h / (avgphiCt * h_t);
        }
    }
}
