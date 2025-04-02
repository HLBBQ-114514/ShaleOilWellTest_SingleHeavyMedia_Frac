using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using ShaleOilWellTest_单重介质压裂封闭边界;
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
        public static double GetPuwD(double u, ReservoirConfig config)
        {
            Debug.WriteLine($"u = {u}, eta = {config.avgeta}, sqrt(u/eta) = {Math.Sqrt(u / config.avgeta)}");
            Func<double, double> f1 = (xD) =>
            {
                return SpecialFunctions.BesselK0(xD);
            };
            Func<double, double> f2 = (xD) =>
            {
                return SpecialFunctions.BesselI0(Sqrt(u / config.avgeta) * xD);
            };
            double intg1 = 2* Sqrt(u / config.avgeta) * GaussLegendreRule.Integrate(f1, 1e-4, Sqrt(u / config.avgeta), 32);//+ GaussLegendreRule.Integrate(f1, 1e-2, 1, 32);
            double intg2 = 2 * GaussLegendreRule.Integrate(f2, 0, 1, 32);

            double omega = config.omega;// config.phi * config.ct * config.h / (config.avgphiCt * config.h_t);//0.5;//
            double zeta = config.zeta;//(config.k / config.mu) * config.h / (config.avgLambda * config.h_t);//1; //

            double value_one = SpecialFunctions.BesselK1(Sqrt(u / config.avgeta) * config.reD);
            double value_two = SpecialFunctions.BesselI1(Sqrt(u / config.avgeta) * config.reD);
            double value_three = 1 / (2 * u * zeta);
            double value_four = config.sf / u;
            double p_D = value_three * value_one / value_two * intg2 + value_three * intg1 + value_four;
            return p_D;

        }
        public static double[] GetPwjD(double u, ReservoirConfig[] config)
        {
            double[] value = new double[2];
            for(int i = 0; i < 2; i++)
            {
                value[i] = 1 / GetPuwD(u, config[i]);
            }
            
            return value;
        }
        public static double GetPwD(double u, ReservoirConfig[] config)
        {
            double pwd = 0;
            foreach(var p in GetPwjD(u, config))
            {
                pwd += p;
            }
            return pwd;
        }

        public static double GetPwCD(double u, ReservoirConfig[] config)
        {
            double value_1 = 0;
            double value_2 = 0;
            for (int i = 0; i < 2; i++)
            {
                value_1 += 1 / GetPuwD(u, config[i]);
            }
            //井筒系数无因次化
            for(int i = 0; i < 2; i++)
            {
                value_2 += config[i].Cs / (6.2832 * config[i].avgphiCt * config[i].h_t * config[i].xf * config[i].xf);
            }
            double value = 1 / (value_1 + u * u * value_2);
            return value;
        }

        public static double[] GetQD(double u, double[] pwd)
        {
            double[] value = new double[2];
            double pw = 0;
            foreach (double p in pwd)
            {
                pw += p;
            }

            for (int i = 0; i < 2; i++)
            {
                value[i] = pw / pwd[i] / u;
            }
            
            return value;
        }

        public static double[] GetQCD(double u, double[] pwd,double pwf)
        {
            double[] value = new double[2];

            for (int i = 0; i < 2; i++)
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
        public static double Getf(double t, int n, ReservoirConfig[] config)
        {
            double f = 0;
            if (n % 2 == 0)
            {
                double ln2 = Math.Log(2) / t;
                Debug.WriteLine("ln2: " + ln2);
                for (int i = 1; i <= n; i++)
                {
                    //f += GetV(n, i) * GetPwD(ln2 * i, config);//无井储                                             
                    f += GetV(n, i) * GetPwCD(ln2 * i, config);//井储
                }

                f = ln2 * f;
            }
            return f;
        }
        public static double[] GetQ(double t, int n, ReservoirConfig[] config)
        {
            double[] qwf = new double[2];
            if (n % 2 == 0)
            {
                double ln2 = Math.Log(2) / t;
                for (int i = 1; i <= n; i++)
                {
                    double[] pwd = GetPwjD(ln2 * i, config);
                    
                    double[] qwd = GetQD(ln2 * i, pwd);
                    for (int j = 0; j < 2; j++)
                    {
                        qwf[j] = GetV(n, i) * qwd[j];
                    }
                }
                for (int i = 0; i < 2; i++)
                {
                    qwf[i] = ln2 * qwf[i];
                }                
            }
            return qwf;
        }

        public static double[] GetCQ(double t, int n, ReservoirConfig[] config)
        {
            double[] qwf = new double[2];
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
                    for (int j = 0; j < 2; j++)
                    {
                        qwf[j] = GetV(n, i) * qwd[j];
                    }
                }
                for (int i = 0; i < 2; i++)
                {
                    qwf[i] = ln2 * qwf[i];
                }
            }
            return qwf;
        }

        public static double GetDf(double pwf, double T, int n,ReservoirConfig[] config)
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
