using MathNet.Numerics;
using MathNet.Numerics.Integration;
using static System.Math;
namespace ShaleOilWellTest
{
    internal class MathUnit
    {
/*        static double omega;
        static double x_j;
        static double s_j;*/

        /*
        * <summary>
        * 在Laplace空间下，求解某个时间的无因次压力(核心)
        * </summary>
        * <param name="u">Laplace变量</param>
        */
        public static double GetPuwD(double u, ReservoirConfigDoubleMedia config)
        {
            //double fu = config.lambdaF * config.lambdaF / (config.lambdaF + u * config.avgLambda*config.omega*config.omegaM/config.kf/config.omegaF)
            //   -config.lambdaF + config.omega * config.omegaF*u/config.zeta;
            double fu = config.lambdaF * config.avgLambda * config.omega * config.omegaM / config.kf / config.omegaF / (config.lambdaF + u * config.avgLambda * config.omega * config.omegaM / config.kf / config.omegaF)
                + config.omega * config.omegaF / config.zeta;
            fu *= u;
            /*Debug.WriteLine("fu:\t"+fu+"\tA:\t"+ config.lambdaF * config.lambdaF / (config.lambdaF + u * config.avgLambda * config.omega * config.omegaM / config.kf / config.omegaF)
                + "\tA2:\t" + u * config.avgLambda * config.omega * config.omegaM / config.kf / config.omegaF
                + "\tB:\t" + config.omega * config.omegaF * u / config.zeta);*/

                //Debug.WriteLine($"u = {u}, eta = {config.avgeta}, sqrt(u/eta) = {Math.Sqrt(u / config.avgeta)}");
            Func<double, double> f1 = (alpha) =>
            {
                return SpecialFunctions.BesselK0(alpha);
            };
            Func<double, double> f2 = (alpha) =>
            {
                return SpecialFunctions.BesselI0(Sqrt(fu)*(Abs(0.732-alpha)));
            };
            double intg1 = 1 / Sqrt(fu) * (GaussLegendreRule.Integrate(f1, 0, Sqrt(fu)*(1+0.732) , 32) + GaussLegendreRule.Integrate(f1, 0, Sqrt(fu)*(1 - 0.732), 32));
            double intg2 = 1 * GaussLegendreRule.Integrate(f2, -1, 1, 32);

            
            double value_one = SpecialFunctions.BesselK1(Sqrt(fu) * config.reD);
            double value_two = SpecialFunctions.BesselI1(Sqrt(fu) * config.reD);
            
            
            //double A = value_one / (u/config.zeta/Sqrt(fu)*(SpecialFunctions.BesselK1(Sqrt(fu)) * value_two-value_one* SpecialFunctions.BesselI1(Sqrt(fu))));
            //double B = A / value_one * value_two;
            double value_three = 1 / (u);
            double value_four = config.sf / u;
            //double p_D = A * intg2 + B * intg1;
            double p_D = value_three * intg1 + value_three * value_one / value_two * intg2 + value_four;
            return p_D;           

        }


        public static double[] GetPwjD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double[] value = new double[config.Length];
            for(int i = 0; i < config.Length; i++)
            {
                value[i] = 1 / GetPuwD(u, config[i]);
            }
            
            return value;
        }

        public static double[] 产量递减GetQ(double u, ReservoirConfigDoubleMedia[] config)
        {
            double[] value = new double[config.Length];
            var pvalue = GetPwjD(u, config);
            for (int i = 0;i < config.Length; i++)
            {
                value[i] = 1 / u / u * pvalue[i];
            }
            return value;
        }
        public static double[] 产量递减GetQLapalce(double t, ReservoirConfigDoubleMedia[] config)
        {
            double[] Q = new double[config.Length];
            if (6 % 2 == 0)
            {
                double ln2 = Math.Log(2) / t;
                //Debug.WriteLine("ln2: " + ln2);
               
                for (int i = 1; i <= 8; i++)
                {
                    var value = 产量递减GetQ(ln2 * i, config);
                    //f += GetV(n, i) * GetPwD(ln2 * i, config);//无井储
                    for(int j = 0; j < Q.Length; j++)
                    {
                        Q[j] += GetV(6, i) * value[j];//井储
                    }         
                }
                for(int i = 0;i < Q.Length; i++)
                {
                    Q[i] *= ln2  * (Math.Log(config[i].reD) - 0.5);
                }
                
            }
            return Q;
        }
        public static double[] 产量递减GetQDi(double[] q,double t)
        {
/*            double[] value = new double[2];
            for (int i = 0; i < q.Rank; i++)
            {
                for(int j = 0; j < q.Length/q.Rank; j++)
                {
                    value[i] += q[i,j];
                }
            }
            return value.Select(q => q / t).ToArray();*/
            return q.Select(qi => qi / t).ToArray();
        }

        public static double[] 产能递减GetdQDi(double[] q, double t, ReservoirConfigDoubleMedia[] config)
        {
            var dq = 产量递减GetQLapalce(t*1.0001, config);
            //dq = 产量递减GetQDi(dq, t);
            for (int i = 0; i < 2; i++)
            {                
                dq[i] = Math.Abs((dq[i] - q[i]) / .0001);
            }
            return dq;
        }

        public static double GetPwD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double pwd = 0;
            foreach(var p in GetPwjD(u, config))
            {
                pwd += 1/p;
            }
            return pwd;
        }
        
        public static double GetPwCD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double value_1 = 0;
            double value_2 = 0;
            for (int i = 0; i < config.Length; i++)
            {
                value_1 += 1 / GetPuwD(u, config[i]);
            }
            //井筒系数无因次化
            for(int i = 0; i < config.Length; i++)
            {
                value_2 += config[i].Cs / (6.2832 * config[i].avgphiCt * config[i].h_t * config[i].xf * config[i].xf);
            }
            double value = 1 / (value_1 + u * u * value_2);
            return value;
        }

        public static double[] GetQD(double u, double[] pwd)
        {
            double[] value = new double[pwd.Length];
            double pw = 0;
            foreach (double p in pwd)
            {
                pw += p;
            }

            for (int i = 0; i < pwd.Length; i++)
            {
                value[i] = pw / pwd[i] / u;
            }
            
            return value;
        }

        public static double[] GetQCD(double u, double[] pwd,double pwf)
        {
            double[] value = new double[pwd.Length];

            for (int i = 0; i < pwd.Length; i++)
            {
                value[i] = pwf * pwd[i] / u;
            }
            return value;
        }
        /*
         * <summary>
         * stefest反演
         * </summary>
         * <param name="t">时间</param>
         * <param name="n">stefest反演系数N</param>
         */
        public static double Getf(double t, int n, ReservoirConfigDoubleMedia[] config)
        {
            double f = 0;
            if (n % 2 == 0)
            {
                double ln2 = Math.Log(2) / t;
                //Debug.WriteLine("ln2: " + ln2);
                for (int i = 1; i <= n; i++)
                {
                    f += GetV(n, i) * GetPwD(ln2 * i, config);//无井储                                             
                    //f += GetV(n, i) * GetPwCD(ln2 * i, config);//井储
                }

                f = ln2 * f;
            }
            return f;
        }
        public static double[] GetQ(double t, int n, ReservoirConfigDoubleMedia[] config)
        {
            double[] qwf = new double[config.Length];
            if (n % 2 == 0)
            {
                double ln2 = Math.Log(2) / t;
                for (int i = 1; i <= n; i++)
                {
                    double[] pwd = GetPwjD(ln2 * i, config);
                    
                    double[] qwd = GetQD(ln2 * i, pwd);
                    for (int j = 0; j < config.Length; j++)
                    {
                        qwf[j] = GetV(n, i) * qwd[j];
                    }
                }
                for (int i = 0; i < config.Length; i++)
                {
                    qwf[i] = ln2 * qwf[i];
                }                
            }
            return qwf;
        }

        public static double[] GetCQ(double t, int n, ReservoirConfigDoubleMedia[] config)
        {
            double[] qwf = new double[config.Length];
            if (n % 2 == 0)
            {
                double ln2 = Math.Log(2) / t;
                for (int i = 1; i <= n; i++)
                {
                    double[] puwd = GetPwjD(ln2 * i, config);
                    //Debug.WriteLine("PwjD: " + puwd[0] + " " + puwd[1]);
                    double pwf = GetPwCD(ln2 * i, config);
                    double[] qwd = GetQCD(ln2 * i, puwd, pwf);
                    //Debug.WriteLine("QD: " + qwd[0] + " " + qwd[1]);
                    for (int j = 0; j < config.Length; j++)
                    {
                        qwf[j] = GetV(n, i) * qwd[j];
                    }
                }
                for (int i = 0; i < config.Length; i++)
                {
                    qwf[i] = ln2 * qwf[i];
                }
            }
            return qwf;
        }

        public static double GetDf(double pwf, double T, int n, ReservoirConfigDoubleMedia[] config)
        {
            double pwf2 = Getf(T * (1.0001), n, config);
            return Math.Abs(pwf2 - pwf) / .0001;

        }

        public static double GetV(int n, int i)
        {
            double V = 0;
            if (n % 2 == 0)
            {
                int k = Convert.ToInt32((i + 1) / 2);
                for (int j = k; j <= Math.Min(n / 2, i); j++)
                {
                    V += GetV_son(n, i, j);
                }
                V = Math.Pow(-1, n / 2 + i) * V;
            }
            return V;
        }

        public static double GetV_son(int n, int i, int k)
        {
            double V_son = 0;
            if (n % 2 == 0)
            {
                V_son = Math.Pow(k, n / 2) * GetFactorial(2 * k) / (GetFactorial(n / 2 - k) * GetFactorial(k) * GetFactorial(k - 1) * GetFactorial(i - k) * GetFactorial(2 * k - i));
            }
            return V_son;

        }

        public static int GetFactorial(int n)
        {
            int N = n;
            if (n >= 2)
            {
                do
                {
                    N *= (n - 1);
                    n--;
                }
                while (n > 2);
            }
            else if (n == 1 || n == 0)
            {
                N = 1;
            }
            else
            {
                N = 0;
            }
            return N;

        }

    }
}
