using MathNet.Numerics;
using MathNet.Numerics.Integration;
using System.Diagnostics;
using static System.Math;
namespace ShaleOilWellTest
{
    internal class MathUnit
    {
         
        #region original welltest code
        /*        static double omega;
                static double x_j;
                static double s_j;*/

        /*
        * <summary>
        * 在Laplace空间下，求解某个时间的无因次压力(核心)
        * </summary>
        * <param name="u">Laplace变量</param>
        */
        public static double GetPuwD_RECTANGULAR(double u, ReservoirConfigDoubleMedia config)
        {
            double xi = config.avgLambda * config.omegaM / config.kf / config.omegaF;
            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * config.avgLambda * config.omega * config.omegaM / config.kf / config.omegaF)
               + config.lambdaF + config.omega * config.omegaF * u / config.zeta;
            /*            
            double fu = config.lambdaF * xi / (config.lambdaF + u * xi)
                    - config.omega * config.omegaF / config.zeta;
            fu *= u;
            */
            
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
            double intg1 = 1 / Sqrt(fu) * (GaussLegendreRule.Integrate(f1, 0, Sqrt(fu)*(config.xfD +0.732) , 32) + GaussLegendreRule.Integrate(f1, 0, Sqrt(fu)*Abs(config.xfD - 0.732), 32));
            double intg2 = GaussLegendreRule.Integrate(f2, -config.xfD, config.xfD, 32);

            
            double value_one = SpecialFunctions.BesselK1(Sqrt(fu) * config.reD);
            double value_two = SpecialFunctions.BesselI1(Sqrt(fu) * config.reD);
            
            
            //double A = value_one / (u/config.zeta/Sqrt(fu)*(SpecialFunctions.BesselK1(Sqrt(fu)) * value_two-value_one* SpecialFunctions.BesselI1(Sqrt(fu))));
            //double B = A / value_one * value_two;
            double value_three = 1 / (u);
            double value_four = config.s / u; ;//
            //double p_D = A * intg2 + B * intg1;
            double p_D = value_three * intg1 + value_three * value_one / value_two * intg2 + value_four;
            return p_D;           

        }

        public static double GetPuwD_CONSTANT(double u, ReservoirConfigDoubleMedia config)
        {

            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * config.avgLambda*config.omega*config.omegaM/config.kf/config.omegaF)
               +config.lambdaF + config.omega * config.omegaF*u/config.zeta;
/*            double xi = config.avgLambda * config.omegaM / config.kf / config.omegaF;
            double fu = config.lambdaF * xi / (config.lambdaF + u * xi)
                        + config.omega * config.omegaF / config.zeta;
            fu *= u;*/
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
                return SpecialFunctions.BesselI0(Sqrt(fu) * (Abs(0.732 - alpha)));
            };
            double intg1 = 1 / Sqrt(fu) * (GaussLegendreRule.Integrate(f1, 0, Sqrt(fu) * (config.xfD + 0.732), 32) + GaussLegendreRule.Integrate(f1, 0, Sqrt(fu) * Abs(config.xfD - 0.732), 32));
            double intg2 = GaussLegendreRule.Integrate(f2, -config.xfD, config.xfD, 32);


            double value_one = SpecialFunctions.BesselK0(Sqrt(fu) * config.reD);
            double value_two = SpecialFunctions.BesselI0(Sqrt(fu) * config.reD);


            //double A = value_one / (u/config.zeta/Sqrt(fu)*(SpecialFunctions.BesselK1(Sqrt(fu)) * value_two-value_one* SpecialFunctions.BesselI1(Sqrt(fu))));
            //double B = A / value_one * value_two;
            double value_three = 1 / (u);
            double value_four = 0;//config.sf / u;
            //double p_D = A * intg2 + B * intg1;
            double p_D = value_three * intg1 - value_three * value_one / value_two * intg2 + value_four;
            return p_D;

        }

        public static double GetPuwD_INFINITY(double u, ReservoirConfigDoubleMedia config)
        {
            double xi = config.avgLambda * config.omegaM * config.omega / config.kf / config.omegaF;
            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * xi)
               + config.lambdaF + config.omega * config.omegaF * u / config.zeta;

            Func<double, double> f1 = (alpha) =>
            {
                return SpecialFunctions.BesselK0(alpha);
            };

            double intg1 = 1 / Sqrt(fu) * (GaussLegendreRule.Integrate(f1, 0, Sqrt(fu) * (config.xfD + 0.732), 32) + GaussLegendreRule.Integrate(f1, 0, Sqrt(fu) * Abs(config.xfD - 0.732), 32));



            //double A = value_one / (u/config.zeta/Sqrt(fu)*(SpecialFunctions.BesselK1(Sqrt(fu)) * value_two-value_one* SpecialFunctions.BesselI1(Sqrt(fu))));
            //double B = A / value_one * value_two;
            double value_three = 1 / (u);
            double value_four = config.s / u;//config.sf / u;
            //double p_D = A * intg2 + B * intg1;
            double p_D =  value_three * intg1 + value_four;
            return p_D;

        }


        public static double GetPwD_cui(double u, ReservoirConfigDoubleMedia config, double r)
        {
            double S = 3;
            double CD = 1e-1;
            var v1 = SpecialFunctions.BesselK0(Sqrt(r * u)) + S * Sqrt(u) * SpecialFunctions.BesselK1(Sqrt(u));
            var demon = u * (Sqrt(u) * SpecialFunctions.BesselK1(Sqrt(u)) + u * CD * (v1));
            return v1 / demon;
        }

        public static double 直线断层边界PwD(double u, ReservoirConfigDoubleMedia config, double r, double op)
        {
            return GetPwD_cui(u, config, 1) + GetPwD_cui(u, config, r);
        }
        #endregion
        #region test code
        public static double TestPwD_Infty(double u, ReservoirConfigDoubleMedia config)
        {
            double xi = config.avgLambda * config.omegaM * config.omega / config.kf / config.omegaF;
            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * xi)  
               + config.lambdaF + config.omega * config.omegaF * u / config.zeta;
            //fu *= u;
            //fu = u;
            int nMax = 100;
            var hD = config.h / config.rw;
            var LD = hD / Math.Cos(config.theta);
            double sinT = Math.Sin(config.theta);
            double cosT = Math.Cos(config.theta);
            Func<double, double> integrand = eta =>
            {
                double rD = Math.Sqrt(eta * eta * Math.Sin(config.theta) * Math.Sin(config.theta));
                double value = SpecialFunctions.BesselK0(
                    rD * Math.Sqrt(fu)
                );
                double sum = 0.0;
                for (int n = 1; n <= nMax; n++)
                {
                    double lambda = Math.Sqrt(
                        fu + n * n * Math.PI * Math.PI / (hD * hD)
                    );

                    sum += 2.0 *
                        SpecialFunctions.BesselK0(rD * lambda) *
                        Math.Cos(n * Math.PI * eta * cosT / hD);
                }
                value += sum;
                return value;
            };

            double integral = GaussLegendreRule.Integrate(
                integrand,
                -LD / 2.0,
                +LD / 2.0,
                64
            );
            return integral / u / LD;

        }

        public static double TestPwD_Rectangular(double u, ReservoirConfigDoubleMedia config)
        {
            double xi = config.avgLambda * config.omegaM * config.omega / config.kf / config.omegaF;
            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * xi)
               + config.lambdaF + config.omega * config.omegaF * u / config.zeta;
            //fu = u;
            int nMax = 100;
            var hD = config.h / config.rw;
            var LD = hD / Math.Cos(config.theta);
            double sinT = Math.Sin(config.theta);
            double cosT = Math.Cos(config.theta);
            Func<double, double> integrand = eta =>
            {
                // 几何距离（与你现有一致）
                double rD = Math.Sqrt(
                    eta * eta * Math.Sin(config.theta) * Math.Sin(config.theta)
                );

                // ---------- n = 0 ----------
                double lambda0 = Math.Sqrt(fu);   // k0 = sqrt(fu)

                double k0_re = lambda0 * config.reD;
                double k0_r = lambda0 * rD;

                double value =
                    SpecialFunctions.BesselK0(k0_r)
                  + (SpecialFunctions.BesselK1(k0_re) / SpecialFunctions.BesselI1(k0_re))
                    * SpecialFunctions.BesselI0(k0_r);

                // ---------- n >= 1 ----------
                double sum = 0.0;

                for (int n = 1; n <= nMax; n++)
                {
                    double lambda = Math.Sqrt(
                        fu + n * n * Math.PI * Math.PI / (hD * hD)
                    );

                    double k_re = lambda * config.reD;
                    double k_r = lambda * rD;

                    double radialKernel =
                        SpecialFunctions.BesselK0(k_r)
                      + (SpecialFunctions.BesselK1(k_re) / SpecialFunctions.BesselI1(k_re))
                        * SpecialFunctions.BesselI0(k_r);

                    sum += 2.0
                        * radialKernel
                        * Math.Cos(n * Math.PI * eta * cosT / hD);
                }

                value += sum;
                return value;
            };

            double integral = GaussLegendreRule.Integrate(
               integrand,
               -LD / 2.0,
               +LD / 2.0,
               64
           );

            return integral / u / LD;
        }

        public static double TestPwD_Constant(double u, ReservoirConfigDoubleMedia config)
        {
            double xi = config.omegaM * config.omega / config.zeta / config.omegaF;
            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * xi)
               + config.lambdaF + config.omega * config.omegaF * u / config.zeta;
            //fu = u;
            int nMax = 100;
            var hD = config.h / config.rw;
            var LD = hD / Math.Cos(config.theta);
            double sinT = Math.Sin(config.theta);
            double cosT = Math.Cos(config.theta);
            Func<double, double> integrand = eta =>
            {
                // 几何距离（与你现有一致）
                double rD = Math.Sqrt(
                    eta * eta * Math.Sin(config.theta) * Math.Sin(config.theta)
                );

                // ---------- n = 0 ----------
                double lambda0 = Math.Sqrt(fu);   // k0 = sqrt(fu)

                double k0_re = lambda0 * config.reD;
                double k0_r = lambda0 * rD;

                double value =
                    SpecialFunctions.BesselK0(k0_r)
                  - (SpecialFunctions.BesselK0(k0_re) / SpecialFunctions.BesselI0(k0_re))
                    * SpecialFunctions.BesselI0(k0_r);

                // ---------- n >= 1 ----------
                double sum = 0.0;

                for (int n = 1; n <= nMax; n++)
                {
                    double lambda = Math.Sqrt(
                       fu + n * n * Math.PI * Math.PI / (hD * hD)
                    );

                    double k_re = lambda * config.reD;
                    double k_r = lambda * rD;

                    double radialKernel =
                        SpecialFunctions.BesselK0(k_r)
                      - (SpecialFunctions.BesselK0(k_re) / SpecialFunctions.BesselI0(k_re))
                        * SpecialFunctions.BesselI0(k_r);

                    sum += 2.0
                        * radialKernel
                        * Math.Cos(n * Math.PI * eta * cosT / hD);
                }

                value += sum;
                return value;
            };

            double integral = GaussLegendreRule.Integrate(
                integrand,
                -LD / 2.0,
                +LD / 2.0,
                32
            );
            return integral / u / LD;

        }
        public static double TestPwCD(double u, ReservoirConfigDoubleMedia config)
        {   
            double PD = TestPwD_Infty(u, config);
            double ct = (config.phiF * config.ctf + config.phiM * config.ctm);
            double CD = config.Cs / (2.0 * Math.PI * ct * config.h_t * config.rw * config.rw);
            CD = 3e-3;
            double S = 1;
            double pwd = (u * PD + S) / (u + CD * u * u * (u * PD + S));
            return pwd;
        }
        #endregion
        #region multilayer test code
        // 按照我其他代码的习惯，把test code 中的代码修改为多层
        public static double MultiTestGetf(double t, int n, ReservoirConfigDoubleMedia[] config)
        {
            return Stehfest.InverseLaplace(t, n, s => MultiTestPwCD(s, config));
        }
        public static double MultiTestGetDf(double pwf, double T, int n, ReservoirConfigDoubleMedia[] config)
        {
            double eps = 1e-3; // 0.1% log step
            double pwf2 = MultiTestGetf(T * (1.0 + eps), n, config);

            return (pwf2 - pwf) / eps;   // 不要 Abs
        }

        public static double MultiTestPwD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double pwc = 0;
            foreach (var c in config)
            {
                pwc += 1/TestPwD_Infty(u, c);
            }
            //pwc /= u * u;
            return 1/pwc;
        }

        public static double MultiTestPwCD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double pwc = 0;
            double CD = 0;
            foreach (var c in config)
            {
                pwc += 1 / TestPwD_Infty(u, c);
                CD += c.Cs / (6.2832 * c.avgphiCt * c.h_t * c.rw * c.rw);
            }
            //pwc /= u * u;            
            
            return 1 / (pwc + u*u*CD);
        }

        #endregion
        #region multilayer SkinFctor test code
        public static double TestPwD_S_Infty(double u, ReservoirConfigDoubleMedia config)
        {
            double xi = config.omegaM * config.omega / config.zeta / config.omegaF;
            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * xi)
               + config.lambdaF + config.omega * config.omegaF * u / config.zeta;
            //fu *= u;
            //fu = u;
            int nMax = 100;
            var hD = config.h / config.rw;
            var LD = hD / Math.Cos(config.theta);
            double sinT = Math.Sin(config.theta);
            double cosT = Math.Cos(config.theta);
            Func<double, double> integrand = eta =>
            {
                // 等效压力点系数（经验值）
                double Cep = 0.732;     // 可作为参数开放
                double rEpD = Cep;     // 已无因次（除以 rw）

                double rD = Math.Sqrt(
                    eta * eta * sinT * sinT + rEpD * rEpD
                );

                double k0 = Math.Sqrt(fu);
                double x0 = k0 * rD;

                double value = SpecialFunctions.BesselK0(x0) 
                + config.s * x0 * SpecialFunctions.BesselK1(x0);

                double sum = 0.0;
                for (int n = 1; n <= nMax; n++)
                {
                    double kn = Math.Sqrt(
                        fu + n * n * Math.PI * Math.PI / (hD * hD)
                    );

                    double x = kn * rD;

                    double term =
                        2.0 *
                        (SpecialFunctions.BesselK0(x) + config.s)
                        * Math.Cos(n * Math.PI * eta * cosT / hD);
                    sum += term;
                }
                value += sum;
                return value;
            };

            double integral = GaussLegendreRule.Integrate(
                integrand,
                -LD / 2.0,
                +LD / 2.0,
                64
            );
            return integral / u / LD;

        }

        public static double TestPwD_S_Rectangular(double u, ReservoirConfigDoubleMedia config)
        {
            double xi = config.avgLambda * config.omegaM * config.omega / config.kf / config.omegaF;
            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * xi)
               + config.lambdaF + config.omega * config.omegaF * u / config.zeta;
            //fu = u;
            int nMax = 100;
            var hD = config.h / config.rw;
            var LD = hD / Math.Cos(config.theta);
            double sinT = Math.Sin(config.theta);
            double cosT = Math.Cos(config.theta);
            Func<double, double> integrand = eta =>
            {
                // 几何距离（与你现有一致）
                double rD = Math.Sqrt(
                    eta * eta * Math.Sin(config.theta) * Math.Sin(config.theta)
                );

                // ---------- n = 0 ----------
                double lambda0 = Math.Sqrt(fu);   // k0 = sqrt(fu)

                double k0_re = lambda0 * config.reD;
                double k0_r = lambda0 * rD;

                double value =
                    SpecialFunctions.BesselK0(k0_r)
                  + (SpecialFunctions.BesselK1(k0_re) / SpecialFunctions.BesselI1(k0_re))
                    * SpecialFunctions.BesselI0(k0_r);

                // ---------- n >= 1 ----------
                double sum = 0.0;

                for (int n = 1; n <= nMax; n++)
                {
                    double lambda = Math.Sqrt(
                        fu + n * n * Math.PI * Math.PI / (hD * hD)
                    );

                    double k_re = lambda * config.reD;
                    double k_r = lambda * rD;

                    double radialKernel =
                        SpecialFunctions.BesselK0(k_r)
                      + (SpecialFunctions.BesselK1(k_re) / SpecialFunctions.BesselI1(k_re))
                        * SpecialFunctions.BesselI0(k_r);

                    sum += 2.0
                        * radialKernel
                        * Math.Cos(n * Math.PI * eta * cosT / hD);
                }

                value += sum;
                return value;
            };

            double integral = GaussLegendreRule.Integrate(
               integrand,
               -LD / 2.0,
               +LD / 2.0,
               64
           );

            return integral / u / LD;
        }

        public static double TestPwD_S_Constant(double u, ReservoirConfigDoubleMedia config)
        {
            double xi = config.omegaM * config.omega / config.zeta / config.omegaF;
            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * xi)
               + config.lambdaF + config.omega * config.omegaF * u / config.zeta;
            //fu = u;
            int nMax = 100;
            var hD = config.h / config.rw;
            var LD = hD / Math.Cos(config.theta);
            double sinT = Math.Sin(config.theta);
            double cosT = Math.Cos(config.theta);
            Func<double, double> integrand = eta =>
            {
                // 几何距离（与你现有一致）
                double rD = Math.Sqrt(
                    eta * eta * Math.Sin(config.theta) * Math.Sin(config.theta)
                );

                // ---------- n = 0 ----------
                double lambda0 = Math.Sqrt(fu);   // k0 = sqrt(fu)

                double k0_re = lambda0 * config.reD;
                double k0_r = lambda0 * rD;

                double value =
                    SpecialFunctions.BesselK0(k0_r)
                  - (SpecialFunctions.BesselK0(k0_re) / SpecialFunctions.BesselI0(k0_re))
                    * SpecialFunctions.BesselI0(k0_r);

                // ---------- n >= 1 ----------
                double sum = 0.0;

                for (int n = 1; n <= nMax; n++)
                {
                    double lambda = Math.Sqrt(
                       fu + n * n * Math.PI * Math.PI / (hD * hD)
                    );

                    double k_re = lambda * config.reD;
                    double k_r = lambda * rD;

                    double radialKernel =
                        SpecialFunctions.BesselK0(k_r)
                      - (SpecialFunctions.BesselK0(k_re) / SpecialFunctions.BesselI0(k_re))
                        * SpecialFunctions.BesselI0(k_r);

                    sum += 2.0
                        * radialKernel
                        * Math.Cos(n * Math.PI * eta * cosT / hD);
                }

                value += sum;
                return value;
            };

            double integral = GaussLegendreRule.Integrate(
                integrand,
                -LD / 2.0,
                +LD / 2.0,
                32
            );
            return integral / u / LD;

        }
        #endregion
        public static double[] GetPwjD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double[] value = new double[config.Length];
            for(int i = 0; i < config.Length; i++)
            {
                value[i] = GetPuwD_RECTANGULAR(u, config[i]);
                //value[i] = GetPuwD_CONSTANT(u, config[i]);
                //value[i] = GetPuwD_INFINITY(u, config[i]);
            }
            
            return value;
        }

        public static double GetPwD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double pwd = 0;
            foreach (var p in GetPwjD(u, config))
            {
                pwd += p;
            }
            return pwd;
        }

        public static double GetPwCD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double value_1 = 0;
            double value_2 = 0;
            for (int i = 0; i < config.Length; i++)
            {
                value_1 += 1 / GetPuwD_RECTANGULAR(u, config[i]);
                //value_1 += 1 / GetPuwD_CONSTANT(u, config[i]);
                //value_1 += 1 / GetPuwD_INFINITY(u, config[i]);;
            }
            //井筒系数无因次化
            for (int i = 0; i < config.Length; i++)
            {
                value_2 += config[i].Cs / (6.2832 * config[i].avgphiCt * config[i].h_t * config[0].xf * config[0].xf);
            }
            double value = 1 / (value_1 + u * u * value_2);
            return value;
        }


        public static double[] 分层产量GetQ(double u, ReservoirConfigDoubleMedia[] config)
        {
            double pwd = 0;
            var pvalue = GetPwjD(u, config).ToArray();
            pwd = GetPwCD(u,config);

            double[] value = new double[config.Length];
            
            for (int i = 0; i < config.Length; i++)
            {
                value[i] = pwd / u / pvalue[i];//* ( config[i].xfD)
            }
            return value;
        }

        public static double[] 产量递减GetQ(double u, ReservoirConfigDoubleMedia[] config)
        {
            double[] value = new double[config.Length];
            var pvalue = GetPwjD(u, config);
            for (int i = 0;i < config.Length; i++)
            {
                value[i] = 1 / u / u / pvalue[i];
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
            return Stehfest.InverseLaplace(t, n, s => GetPwCD(s, config));
            /*double f = 0;
            if (n % 2 == 0)
            {
                double ln2 = Math.Log(2) / t;
                //Debug.WriteLine("ln2: " + ln2);
                for (int i = 1; i <= n; i++)
                {
                    //f += GetV(n, i) * GetPwD(ln2 * i, config);//无井储                                             
                      f += GetV(n, i) * GetPwCD(ln2 * i, config);//井储
                }

                f = ln2 * f;
            }
            return f;*/
        }

        public static double Getf_cui(double t, int n, ReservoirConfigDoubleMedia config)
        {   
            return Stehfest.InverseLaplace(t, n, s => TestPwD_S_Infty(s, config));
            //return Stehfest.InverseLaplace(t, n, s => TestPwD_Constant(s, config));
            //return Stehfest.InverseLaplace(t, n, s => TestPwD_Rectangular(s, config));
            //return Stehfest.InverseLaplace(t, n, s => 直线断层边界PwD(s, config, 10000,1));
        }
        public static double GetDf_cui(double pwf, double T, int n, ReservoirConfigDoubleMedia config)
        {
/*            double pwf2 = Getf_cui(T * (1.0001), n, config);
            return Math.Abs(pwf2 - pwf) / .0001;*/
            double eps = 1e-3; // 0.1% log step
            double pwf2 = Getf_cui(T * (1.0 + eps), n, config);

            return (pwf2 - pwf) / eps;   // 不要 Abs
        }
        public static double[] GetQ(double t, int n, ReservoirConfigDoubleMedia[] config)
        {
            double[] qwf = new double[config.Length];
            if (n % 2 == 0)
            {
                double ln2 = Math.Log(2) / t;
                for (int i = 1; i <= n; i++)
                {
                    double[] qwd = 分层产量GetQ(ln2 * i, config);                   
        
                    for (int j = 0; j < config.Length; j++)
                    {
                        qwf[j] += GetV(n, i) * qwd[j];
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
            double pwf2 = Getf(T * (1.001), n, config);
            return Math.Abs(pwf2 - pwf) / .001;

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
