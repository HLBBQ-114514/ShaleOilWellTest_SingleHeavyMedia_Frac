using MathNet.Numerics;
using MathNet.Numerics.Integration;
using System;
using static System.Math;
namespace ShaleOilWellTest
{
    internal class MathUnit
    {
        private static class BesselSafe
        {
            private const double Euler = 0.5772156649015328606;
            private const double Small = 1e-6;
            private const double Large = 50.0;

            public static double K0(double x)
            {
                if (x < Small)
                {
                    double y = x * x / 4.0;
                    double logTerm = -Log(x / 2.0) - Euler;
                    return logTerm + y * (1 - logTerm) / 2.0;
                }

                if (x > Large)
                {
                    double factor = Sqrt(PI / (2.0 * x));
                    return factor * Exp(-x) * (1.0 + 1.0 / (8.0 * x));
                }

                var v = SpecialFunctions.BesselK0(x);
                if (double.IsFinite(v)) return v;
                double fallback = Sqrt(PI / (2.0 * x)) * Exp(-x);
                return fallback;
            }

            public static double K1(double x)
            {
                if (x < Small)
                {
                    // K1 ~ 1/x + x/2*(ln(x/2)+gamma-1/2)
                    double logTerm = Log(x / 2.0) + Euler - 0.5;
                    return 1.0 / x + x * logTerm / 2.0;
                }

                if (x > Large)
                {
                    double factor = Sqrt(PI / (2.0 * x)) * Exp(-x);
                    return factor * (1.0 + 3.0 / (8.0 * x));
                }

                var v = SpecialFunctions.BesselK1(x);
                if (double.IsFinite(v)) return v;
                double fallback = Sqrt(PI / (2.0 * x)) * Exp(-x);
                return fallback;
            }

            public static double I0(double x)
            {
                double ax = Abs(x);
                if (ax < Small)
                {
                    double y = x * x / 4.0;
                    return 1.0 + y + y * y / 4.0;
                }

                if (ax > Large)
                {
                    double factor = Exp(ax) / Sqrt(2.0 * PI * ax);
                    return factor * (1.0 + 1.0 / (8.0 * ax));
                }

                var v = SpecialFunctions.BesselI0(x);
                if (double.IsFinite(v)) return v;
                double fallback = Exp(ax) / Sqrt(2.0 * PI * ax);
                return fallback;
            }

            public static double I1(double x)
            {
                double ax = Abs(x);
                if (ax < Small)
                {
                    return x / 2.0 + x * x * x / 16.0;
                }

                if (ax > Large)
                {
                    double factor = Exp(ax) / Sqrt(2.0 * PI * ax);
                    double series = 1.0 - 3.0 / (8.0 * ax);
                    double val = factor * series;
                    return x < 0 ? -val : val;
                }

                var v = SpecialFunctions.BesselI1(x);
                if (double.IsFinite(v)) return v;
                double fallback = Exp(ax) / Sqrt(2.0 * PI * ax);
                return x < 0 ? -fallback : fallback;
            }
        }

        // 数值稳定的 Bessel 评估：大 x 用渐近，极小 x 用级数，不做“冻结”裁剪
        private static double K0Scaled(double x)
        {
            if (x <= 0) return double.NaN;
            if (x < 1e-4)
            {
                double y = x * x / 4.0;
                double logTerm = -Log(x / 2.0) - 0.5772156649015328606;
                return logTerm + y * (1 - logTerm) / 2.0;
            }
            if (x > 80)
            {
                double pref = Sqrt(PI / (2.0 * x));
                return pref * Exp(-x) * (1.0 + 1.0 / (8.0 * x));
            }
            return BesselSafe.K0(x);
        }

        private static double K1Scaled(double x)
        {
            if (x <= 0) return double.NaN;
            if (x < 1e-4)
            {
                double logTerm = Log(x / 2.0) + 0.5772156649015328606 - 0.5;
                return 1.0 / x + x * logTerm / 2.0;
            }
            if (x > 80)
            {
                double pref = Sqrt(PI / (2.0 * x));
                return pref * Exp(-x) * (1.0 + 3.0 / (8.0 * x));
            }
            return BesselSafe.K1(x);
        }

        private static double I0Scaled(double x)
        {
            double ax = Abs(x);
            if (ax < 1e-4)
            {
                double y = x * x / 4.0;
                return 1.0 + y + y * y / 4.0;
            }
            if (ax > 80)
            {
                double pref = Exp(ax) / Sqrt(2.0 * PI * ax);
                return pref * (1.0 + 1.0 / (8.0 * ax));
            }
            return BesselSafe.I0(x);
        }

        private static double I1Scaled(double x)
        {
            double ax = Abs(x);
            if (ax < 1e-4)
            {
                return x / 2.0 + x * x * x / 16.0;
            }
            if (ax > 80)
            {
                double pref = Exp(ax) / Sqrt(2.0 * PI * ax);
                double series = 1.0 - 3.0 / (8.0 * ax);
                double val = pref * series;
                return x < 0 ? -val : val;
            }
            return BesselSafe.I1(x);
        }

        private static double LogI0(double x)
        {
            double ax = Abs(x);
            if (ax > 50)
            {
                // I0 ~ exp(x)/sqrt(2πx)
                return ax - 0.5 * Log(2 * PI * ax);
            }
            if (ax < 1e-4)
            {
                double y = x * x / 4.0; // series 1 + y + y^2/4
                double val = 1.0 + y + y * y / 4.0;
                return Log(val);
            }
            return Log(BesselSafe.I0(x));
        }

        private static double LogI1(double x)
        {
            double ax = Abs(x);
            if (ax > 50)
            {
                // I1 ~ exp(x)/sqrt(2πx) * (1 - 3/(8x))
                double leading = ax - 0.5 * Log(2 * PI * ax);
                double corr = Log(1.0 - 3.0 / (8.0 * ax));
                double logVal = leading + corr;
                return x < 0 ? logVal + Log(-1.0) : logVal; // sign handled by caller
            }
            if (ax < 1e-4)
            {
                double val = x / 2.0 + x * x * x / 16.0;
                return Log(Abs(val));
            }
            return Log(Abs(BesselSafe.I1(x)));
        }

        private static double LogK0(double x)
        {
            if (x < 1e-4)
            {
                double logTerm = -Log(x / 2.0) - 0.5772156649015328606;
                // K0 ≈ logTerm
                return Log(logTerm);
            }
            if (x > 50)
            {
                // K0 ~ sqrt(pi/(2x)) * exp(-x)
                return -x + 0.5 * (Log(PI) - Log(2 * x));
            }
            return Log(BesselSafe.K0(x));
        }

        private static double LogK1(double x)
        {
            if (x < 1e-4)
            {
                // K1 ~ 1/x
                return -Log(x);
            }
            if (x > 50)
            {
                // K1 ~ sqrt(pi/(2x)) * exp(-x) * (1 + 3/(8x))
                double leading = -x + 0.5 * (Log(PI) - Log(2 * x));
                double corr = Log(1.0 + 3.0 / (8.0 * x));
                return leading + corr;
            }
            return Log(BesselSafe.K1(x));
        }

        private static double RatioK1OverI1(double x)
        {
            if (x <= 0) return double.NaN;
            double logK = LogK1(x);
            double logI = LogI1(x);
            return Exp(logK - logI);
        }

        private static double RatioK0OverI0(double x)
        {
            if (x <= 0) return double.NaN;
            double logK = LogK0(x);
            double logI = LogI0(x);
            return Exp(logK - logI);
        }

        /// <summary>
        /// 通用 fu 计算，供不同边界条件复用。
        /// </summary>
        private static double ComputeFu(double u, ReservoirConfigDoubleMedia config)
        {
            double denominator = config.lambdaF + u * config.avgLambda * config.omega * config.omegaM / config.kf / config.omegaF;
            double fuRaw = -config.lambdaF * config.lambdaF / denominator
                   + config.lambdaF
                   + config.omega * config.omegaF * u / config.zeta;
            return fuRaw;
        }

        /// <summary>
        /// 计算裂缝两项积分：intg1（K0）与 intg2（I0）。
        /// </summary>
        private static (double intg1, double intg2) ComputeFractureIntegrals(double sqrtFu, ReservoirConfigDoubleMedia config)
        {
            Func<double, double> f1 = alpha => K0Scaled(alpha);
            Func<double, double> f2 = alpha => I0Scaled(sqrtFu * Abs(0.732 - alpha));

            double intg1 = 1 / sqrtFu * (
                GaussLegendreRule.Integrate(f1, 0, sqrtFu * (config.xfD + 0.732), 32) +
                GaussLegendreRule.Integrate(f1, 0, sqrtFu * Abs(config.xfD - 0.732), 32));

            double intg2 = GaussLegendreRule.Integrate(f2, -config.xfD, config.xfD, 32);
            return (intg1, intg2);
        }


        #region original welltest code

        /// <summary>
        /// 在Laplace空间下，求解某个时间的无因次压力(核心)
        /// </summary>
        /// <param name="u">Laplace变量</param>
        /// <param name="config"></param>
        /// <returns></returns>
        public static double GetPuwD_RECTANGULAR(double u, ReservoirConfigDoubleMedia config)
        {
            double fu = ComputeFu(u, config);
            if (fu <= 0) return double.NaN;
            double sqrtFu = Sqrt(fu);
            //Debug.WriteLine($"u = {u}, eta = {config.avgeta}, sqrt(u/eta) = {Math.Sqrt(u / config.avgeta)}");
            var (intg1, intg2) = ComputeFractureIntegrals(sqrtFu, config);

            double value_three = 1 / u;
            double value_four = config.s * value_three; ;//
            double ratio = RatioK1OverI1(sqrtFu * config.reD);
            double p_D = value_three * intg1 + value_three * ratio * intg2 + value_four;
            return p_D;

        }

        public static double GetPuwD_CONSTANT(double u, ReservoirConfigDoubleMedia config)
        {

            double fu = ComputeFu(u, config);
            if (fu <= 0) return double.NaN;
            double sqrtFu = Sqrt(fu);
            //Debug.WriteLine($"u = {u}, eta = {config.avgeta}, sqrt(u/eta) = {Math.Sqrt(u / config.avgeta)}");
            var (intg1, intg2) = ComputeFractureIntegrals(sqrtFu, config);

            double value_three = 1 / u;
            double value_four = 0;//config.sf / u;
            double ratio = RatioK0OverI0(sqrtFu * config.reD);
            double p_D = value_three * intg1 - value_three * ratio * intg2 + value_four;
            return p_D;

        }

        public static double GetPuwD_INFINITY(double u, ReservoirConfigDoubleMedia config)
        {
            double fu = ComputeFu(u, config);
            if (fu <= 0) return double.NaN;
            double sqrtFu = Sqrt(fu);

            var (intg1, _) = ComputeFractureIntegrals(sqrtFu, config);

            double value_three = 1 / u;
            double value_four = config.s * value_three;//config.sf / u;

            double p_D = value_three * intg1 + value_four;
            return p_D;

        }


        public static double GetPwD_cui(double u, ReservoirConfigDoubleMedia config, double r)
        {
            double S = 3;
            double CD = 1e-1;
            double sqrtRU = Sqrt(r * u);
            double sqrtU = Sqrt(u);
            var v1 = K0Scaled(sqrtRU) + S * sqrtU * K1Scaled(sqrtU);
            var demon = u * (sqrtU * K1Scaled(sqrtU) + u * CD * v1);
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
            if (fu <= 0) return double.NaN;
            int nMax = 100;
            var hD = config.h / config.rw;
            var LD = hD / Math.Cos(config.theta);
            double sinT = Math.Sin(config.theta);
            double cosT = Math.Cos(config.theta);
            Func<double, double> integrand = eta =>
            {
                double rD = Math.Sqrt(eta * eta * Math.Sin(config.theta) * Math.Sin(config.theta));
                double value = K0Scaled(rD * Math.Sqrt(fu));
                double sum = 0.0;
                for (int n = 1; n <= nMax; n++)
                {
                    double lambda = Math.Sqrt(
                        fu + n * n * Math.PI * Math.PI / (hD * hD)
                    );

                    sum += 2.0 *
                        K0Scaled(rD * lambda) *
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
            if (fu <= 0) return double.NaN;
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
                    K0Scaled(k0_r)
                  + RatioK1OverI1(k0_re)
                    * I0Scaled(k0_r);

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
                        K0Scaled(k_r)
                      + RatioK1OverI1(k_re)
                        * I0Scaled(k_r);

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
            if (fu <= 0) return double.NaN;
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
                    K0Scaled(k0_r)
                  - RatioK0OverI0(k0_re)
                    * I0Scaled(k0_r);

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
                        K0Scaled(k_r)
                      - RatioK0OverI0(k_re)
                        * I0Scaled(k_r);

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
            double PD = TestPwD_S_Constant(u, config);
            double ct = (config.phiF * config.ctf + config.phiM * config.ctm);
            double CD = config.Cs / (2.0 * Math.PI * ct * config.h_t * config.rw * config.rw);
            CD = 5e-3;
            double S = 0;
            double pwd = (u * PD + S) / (u + CD * u * u * (u * PD + S));
            return pwd;
        }
        #endregion
        #region multilayer test code
        // 按照我其他代码的习惯，把test code 中的代码修改为多层
        public static double MultiTestGetf(double t, int n, ReservoirConfigDoubleMedia[] config)
        {
            return Stehfest.InverseLaplace(t, n, s => MultiTestPwCD_Original(s, config));
            //return Stehfest.InverseLaplace(t, n, s => MultiTestPwCD(s, config));
        }
        public static double MultiTestGetDf(double pwf, double T, int n, ReservoirConfigDoubleMedia[] config)
        {
            double eps = Max(1e-5, Min(5e-3, 0.01 / (1.0 + T)));
            double tPlus = T * (1.0 + eps);
            double tMinus = T * (1.0 - eps);
            double fPlus = MultiTestGetf(tPlus, n, config);
            double fMinus = MultiTestGetf(tMinus, n, config);

            return (fPlus - fMinus) / (2.0 * eps);
        }

        public static double MultiTestPwD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double numerator = 1.0 / u;
            double denominator = 0.0;

            foreach (var c in config)
            {
                double G = TestPwD_S_Constant(u, c); // 已积分的 \bar{G}_{wDj}(u)
                if (double.IsNaN(G) || double.IsInfinity(G)) continue;

                double pi = c.PiD;
                double denomTerm = (u * G - pi);
                if (denomTerm == 0) continue;

                denominator += 1.0 / denomTerm;

                double numTerm = (u * u * G - pi);
                if (numTerm == 0) continue;
                numerator += pi / numTerm;
            }

            if (denominator == 0) return double.NaN;
            return numerator / denominator;
        }
        public static double MultiTestPwCD_Original(double u, ReservoirConfigDoubleMedia[] config)
        {
            double pwc = 0;
            double CD = 0;
            foreach (var c in config)
            {
                pwc += 1 / TestPwD_S_Constant_MathNet(u, c);//TestPwD_S_Constant_MathNet
                CD += c.Cs / (6.2832 * c.avgphiCt * c.h_t * c.rw * c.rw);
            }

            return 1 / (pwc + u * u * CD);
        }

        public static double MultiTestPwCD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double sumCD = 0.0;
            double numerator = 1.0 / u;
            double denominator = 0.0;
            double[] piD = new double[config.Length];
            foreach (var c in config)
            {
                double G = TestPwD_S_Constant(u, c); // \bar{G}_{wDj}(u)
                if (double.IsNaN(G) || double.IsInfinity(G)) continue;

                double pi = 0;
                double CDj = c.Cs / (6.2832 * c.avgphiCt * c.h_t * c.rw * c.rw);
                sumCD += CDj;

                double denomTerm = (u * G - pi);
                if (denomTerm != 0)
                {
                    denominator += 1.0 / denomTerm;
                }

                double numTerm = (u * u * G - pi);
                double invNum = numTerm != 0 ? 1.0 / numTerm : 0.0;
                numerator += pi * (invNum + CDj);
            }

            denominator += u * sumCD;
            if (denominator == 0) return double.NaN;
            return numerator / denominator;
        }

        #endregion
        #region multilayer SkinFctor test code
        public static double TestPwD_S_Infty(double u, ReservoirConfigDoubleMedia config)
        {
            double xi = config.omegaM * config.omega / config.zeta / config.omegaF;
            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * xi)
               + config.lambdaF + config.omega * config.omegaF * u / config.zeta;
            if (fu <= 0) return double.NaN;
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

                double value = K0Scaled(x0) 
                + config.s * x0 * K1Scaled(x0);

                double sum = 0.0;
                for (int n = 1; n <= nMax; n++)
                {
                    double kn = Math.Sqrt(
                        fu + n * n * Math.PI * Math.PI / (hD * hD)
                    );

                    double x = kn * rD;

                    double term =
                        2.0 *
                        (K0Scaled(x) + config.s * x * K1Scaled(x))
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
            if (fu <= 0) return double.NaN;
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
                    K0Scaled(k0_r)
                  + RatioK1OverI1(k0_re) * I0Scaled(k0_r)
                  + config.s * k0_r *
                    (K1Scaled(k0_r)
                     - RatioK1OverI1(k0_re) * I1Scaled(k0_r));

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
                        K0Scaled(k_r)
                      + RatioK1OverI1(k_re) * I0Scaled(k_r)
                      + config.s * k_r *
                        (K1Scaled(k_r)
                         - RatioK1OverI1(k_re) * I1Scaled(k_r));

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
            if (fu <= 0) return double.NaN;
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
                    K0Scaled(k0_r)
                  - RatioK0OverI0(k0_re) * I0Scaled(k0_r)
                  + config.s * k0_r *
                    (K1Scaled(k0_r)
                     - RatioK0OverI0(k0_re) * I1Scaled(k0_r));

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
                        K0Scaled(k_r)
                      - RatioK0OverI0(k_re) * I0Scaled(k_r)
                      + config.s * k_r *
                        (K1Scaled(k_r)
                         - RatioK0OverI0(k_re) * I1Scaled(k_r));

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

        /// <summary>
        /// 原始 MathNet Bessel 版本（无缩放），用于对比数值差异。
        /// </summary>
        public static double TestPwD_S_Constant_MathNet(double u, ReservoirConfigDoubleMedia config)
        {
            double xi = config.omegaM * config.omega / config.zeta / config.omegaF;
            double fu = -config.lambdaF * config.lambdaF / (config.lambdaF + u * xi)
               + config.lambdaF + config.omega * config.omegaF * u / config.zeta;
            if (fu <= 0) return double.NaN;
            int nMax = 100;
            var hD = config.h / config.rw;
            var LD = hD / Math.Cos(config.theta);
            double sinT = Math.Sin(config.theta);
            double cosT = Math.Cos(config.theta);

            Func<double, double> integrand = eta =>
            {
                double rD = Math.Sqrt(
                    eta * eta * sinT * sinT
                );

                double lambda0 = Math.Sqrt(fu);
                double k0_re = lambda0 * config.reD;
                double k0_r = lambda0 * rD;

                double value =
                    SpecialFunctions.BesselK0(k0_r)
                  - (SpecialFunctions.BesselK0(k0_re) / SpecialFunctions.BesselI0(k0_re))
                    * SpecialFunctions.BesselI0(k0_r)
                  + config.s * k0_r *
                    (SpecialFunctions.BesselK1(k0_r)
                     - (SpecialFunctions.BesselK0(k0_re) / SpecialFunctions.BesselI0(k0_re))
                        * SpecialFunctions.BesselI1(k0_r));

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
                        * SpecialFunctions.BesselI0(k_r)
                      + config.s * k_r *
                        (SpecialFunctions.BesselK1(k_r)
                         - (SpecialFunctions.BesselK0(k_re) / SpecialFunctions.BesselI0(k_re))
                            * SpecialFunctions.BesselI1(k_r));

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
        #endregion
        public static double[] GetPwjD(double u, ReservoirConfigDoubleMedia[] config)
        {
            double[] value = new double[config.Length];
            for (int i = 0; i < config.Length; i++)
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

        public static double GetPwCD(double u, ReservoirConfigDoubleMedia[] config, double[]? precomputedPuwD = null)
        {
            double value_1 = 0;
            double value_2 = 0;
            if (precomputedPuwD != null)
            {
                if (precomputedPuwD.Length != config.Length)
                {
                    throw new ArgumentException("预计算的压力长度必须与层数一致", nameof(precomputedPuwD));
                }

                for (int i = 0; i < precomputedPuwD.Length; i++)
                {
                    value_1 += 1 / precomputedPuwD[i];
                }
            }
            else
            {
                for (int i = 0; i < config.Length; i++)
                {
                    value_1 += 1 / GetPuwD_RECTANGULAR(u, config[i]);
                    //value_1 += 1 / GetPuwD_CONSTANT(u, config[i]);
                    //value_1 += 1 / GetPuwD_INFINITY(u, config[i]); ;
                }
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
            pwd = GetPwCD(u, config, pvalue);

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
            for (int i = 0; i < config.Length; i++)
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
                    for (int j = 0; j < Q.Length; j++)
                    {
                        //Q[j] += GetV(6, i) * value[j];//井储
                    }
                }
                for (int i = 0; i < Q.Length; i++)
                {
                    Q[i] *= ln2 * (Math.Log(config[i].reD) - 0.5);
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

        public static double[] GetQCD(double u, double[] pwd, double pwf)
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
            return Stehfest.InverseLaplace(t, n, s => TestPwCD(s, config));//TestPwCD
            //return Stehfest.InverseLaplace(t, n, s => TestPwD_Constant(s, config));
            //return Stehfest.InverseLaplace(t, n, s => TestPwD_Rectangular(s, config));
            //return Stehfest.InverseLaplace(t, n, s => 直线断层边界PwD(s, config, 10000,1));
        }
        public static double GetDf_cui(double pwf, double T, int n, ReservoirConfigDoubleMedia config)
        {
            // 对数刻度对称差分，进一步压制导数噪声
            double eps = Max(1e-5, Min(5e-3, 0.01 / (1.0 + T)));
            double tPlus = T * (1.0 + eps);
            double tMinus = T * (1.0 - eps);
            double fPlus = Getf_cui(tPlus, n, config);
            double fMinus = Getf_cui(tMinus, n, config);

            return (fPlus - fMinus) / (2.0 * eps);
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
                        //qwf[j] += GetV(n, i) * qwd[j];
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
                    double pwf = GetPwCD(ln2 * i, config, puwd);
                    double[] qwd = GetQCD(ln2 * i, puwd, pwf);
                    //Debug.WriteLine("QD: " + qwd[0] + " " + qwd[1]);
                    for (int j = 0; j < config.Length; j++)
                    {
                        //qwf[j] = GetV(n, i) * qwd[j];
                    }
                }
                for (int i = 0; i < config.Length; i++)
                {
                    qwf[i] = ln2 * qwf[i];
                }
            }
            return qwf;
        }
    }
}
