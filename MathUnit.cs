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
                return SpecialFunctions.BesselI0( xD);
            };
            double intg1 = 2/ Sqrt(u / config.avgeta) * GaussLegendreRule.Integrate(f1, 1e-4, Sqrt(u / config.avgeta), 32);//+ GaussLegendreRule.Integrate(f1, 1e-2, 1, 32);
            double intg2 = 2 * Sqrt(u / config.avgeta) * GaussLegendreRule.Integrate(f2, 1e-4, Sqrt(u / config.avgeta), 32);

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

        /// <summary>
        /// 变形第一类整数阶贝塞尔函数I
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double TransformativeIntegerBessel1stFunction(int n, double x)
        {
            int i, m;
            double t, y, p = 0, b0, b1, q;
            double[] a ={ 1.0,3.5156229,3.0899424,1.2067492,
                                 0.2659732,0.0360768,0.0045813};
            double[] b ={ 0.5,0.87890594,0.51498869,
                      0.15084934,0.02658773,0.00301532,0.00032411};
            double[] c ={ 0.39894228,0.01328592,0.00225319,
                                -0.00157565,0.00916281,-0.02057706,
                                 0.02635537,-0.01647633,0.00392377};
            double[] d ={ 0.39894228,-0.03988024,-0.00362018,
                                0.00163801,-0.01031555,0.02282967,
                                -0.02895312,0.01787654,-0.00420059};
            if (n < 0) n = -n;
            t = Math.Abs(x);
            if (n != 1)
            {
                if (t < 3.75)
                {
                    y = (x / 3.75) * (x / 3.75);
                    p = a[6];
                    for (i = 5; i >= 0; i--) p = p * y + a[i];
                }
                else
                {
                    y = 3.75 / t;
                    p = c[8];
                    for (i = 7; i >= 0; i--) p = p * y + c[i];
                    p = p * Math.Exp(t) / Math.Sqrt(t);
                }
            }
            if (n == 0) return (p);
            q = p;
            if (t < 3.75)
            {
                y = (x / 3.75) * (x / 3.75);
                p = b[6];
                for (i = 5; i >= 0; i--) p = p * y + b[i];
                p = p * t;
            }
            else
            {
                y = 3.75 / t;
                p = d[8];
                for (i = 7; i >= 0; i--) p = p * y + d[i];
                p = p * Math.Exp(t) / Math.Sqrt(t);
            }
            if (x < 0.0) p = -p;
            if (n == 1) return (p);
            if (x == 0.0) return (0.0);
            y = 2.0 / t;
            t = b0 = 0.0;
            b1 = 1.0;
            m = n + (int)Math.Sqrt(40.0 * n);
            m = 2 * m;
            for (i = m; i > 0; i--)
            {
                p = b0 + i * y * b1;
                b0 = b1;
                b1 = p;
                if (Math.Abs(b1) > 1.0e+10)
                {
                    t = t * 1.0e-10;
                    b0 = b0 * 1.0e-10;
                    b1 = b1 * 1.0e-10;
                }
                if (i == n) t = b0;
            }
            p = t * q / b1;
            if ((x < 0.0) && (n % 2 == 1)) p = -p;
            return (p);
        }

        /// <summary>
        /// 变形第二类整数阶贝塞尔函数K
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double TransformativeIntegerBessel2ndFunction(int n, double x)
        {
            int i;
            double y, p = 0, b0, b1;

            double[] a ={ -0.57721566,0.4227842,0.23069756,
           0.0348859,0.00262698,0.0001075,0.0000074};
            double[] b ={ 1.0,0.15443144,-0.67278579,
           -0.18156897,-0.01919402,-0.00110404,-0.00004686};
            double[] c ={ 1.25331414,-0.07832358,0.02189568,
           -0.01062446,0.00587872,-0.0025154,0.00053208};
            double[] d ={ 1.25331414,0.23498619,-0.0365562,
           0.01504268,-0.00780353,0.00325614,-0.00068245};
            if (n < 0) n = -n;
            if (x < 0.0) x = -x;
            if (x == 0.0) return (1.0e+70);
            if (n != 1)
            {
                if (x < 2.0 || x == 2.0)
                {
                    y = x * x / 4.0; p = a[6];
                    for (i = 5; i >= 0; i--) p = p * y + a[i];
                    p = p - TransformativeIntegerBessel1stFunction(0, x) * Math.Log(x / 2.0);
                }
                else
                {
                    y = 2.0 / x; p = c[6];
                    for (i = 5; i >= 0; i--) p = p * y + c[i];
                    p = p * Math.Exp(-x) / Math.Sqrt(x);
                }
            }
            if (n == 0) return (p);
            b0 = p;
            if (x < 2.0 || x == 2.0)
            {
                y = x * x / 4.0;
                p = b[6];
                for (i = 5; i >= 0; i--) p = p * y + b[i];
                p = p / x + TransformativeIntegerBessel1stFunction(1, x) * Math.Log(x / 2.0);
            }
            else
            {
                y = 2.0 / x;
                p = d[6];
                for (i = 5; i >= 0; i--) p = p * y + d[i];
                p = p * Math.Exp(-x) / Math.Sqrt(x);
            }
            if (n == 1) return (p);
            b1 = p;
            y = 2.0 / x;
            for (i = 1; i < n; i++)
            {
                p = b0 + i * y * b1;
                b0 = b1;
                b1 = p;
            }
            return (p);
        }

    }
}

