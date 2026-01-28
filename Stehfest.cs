using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Concurrent;

public static class Stehfest
{
    // ===== 1) log 阶乘表（够用到 n<=64 一般都稳）=====
    // Stehfest 需要用到 LogFactorial(2k) 最高到 2*(n/2)=n
    // 但中间项会有 2k, 2k-i 等，建议预留更大一点：4*n
    private static readonly double[] logFact = InitLogFact(2000);

    private static double[] InitLogFact(int maxN)
    {
        double[] f = new double[maxN + 1];
        f[0] = 0.0;
        for (int i = 1; i <= maxN; i++) f[i] = f[i - 1] + Math.Log(i);
        return f;
    }

    private static double LogFactorial(int n)
    {
        if (n <= 1) return 0.0;
        if (n < logFact.Length) return logFact[n];

        // Stirling (足够用于超出预计算范围)
        double x = n;
        return x * Math.Log(x) - x + 0.5 * Math.Log(2.0 * Math.PI * x)
               + 1.0 / (12.0 * x) - 1.0 / (360.0 * x * x * x);
    }

    // ===== 2) V 系数缓存：n -> V[1..n]（下标从 1 用起来更直观）=====
    private static readonly ConcurrentDictionary<int, double[]> VCache = new();

    /// <summary>
    /// 获取 n 阶 Stehfest 系数数组 V（长度 n+1，下标 1..n 有效）
    /// </summary>
    public static double[] GetVCoefficients(int n)
    {
        if (n <= 0 || (n % 2) != 0) throw new ArgumentException("n 必须为正偶数");
        return VCache.GetOrAdd(n, ComputeV);
    }

    // ===== 3) 用 log-sum-exp 计算每个 V_i（最关键的稳定性改进）=====
    private static double[] ComputeV(int n)
    {
        int half = n / 2;
        var V = new double[n + 1]; // V[0] unused

        // 预先计算 log(k) 以减少重复
        double[] logK = new double[half + 1];
        for (int k = 1; k <= half; k++) logK[k] = Math.Log(k);

        for (int i = 1; i <= n; i++)
        {
            int kMin = (i + 1) / 2;
            int kMax = Math.Min(half, i);

            // log-sum-exp：sum_k exp(logV_k)
            double maxLog = double.NegativeInfinity;
            double[] logs = new double[kMax - kMin + 1];

            int idx = 0;
            for (int k = kMin; k <= kMax; k++)
            {
                // logV = (n/2)*ln(k) + ln((2k)!) - ln((n/2-k)!) - ln(k!) - ln((k-1)!)
                //        - ln((i-k)!) - ln((2k-i)!)
                double logV =
                    half * logK[k]
                    + LogFactorial(2 * k)
                    - LogFactorial(half - k)
                    - LogFactorial(k)
                    - LogFactorial(k - 1)
                    - LogFactorial(i - k)
                    - LogFactorial(2 * k - i);

                logs[idx++] = logV;
                if (logV > maxLog) maxLog = logV;
            }

            // 计算 sumExp = Σ exp(logV - maxLog)
            double sumExp = 0.0;
            for (int j = 0; j < logs.Length; j++)
                sumExp += Math.Exp(logs[j] - maxLog);

            double sum = Math.Exp(maxLog) * sumExp;

            // 符号项：(-1)^{n/2 + i}
            double sign = (((half + i) & 1) == 0) ? 1.0 : -1.0;
            V[i] = sign * sum;
        }

        return V;
    }

    // ===== 4) 标量反演（高性能：系数一次取出，循环内不再算 V）=====
    public static double InverseLaplace(double t, int n, Func<double, double> F, bool skipNaNInf = true)
    {
        if (t <= 0) throw new ArgumentException("t 必须 > 0");
        var V = GetVCoefficients(n);

        double ln2 = Math.Log(2.0);
        double factor = ln2 / t;

        // Kahan summation
        double sum = 0.0, c = 0.0;

        for (int i = 1; i <= n; i++)
        {
            double s = i * factor;
            double fs = F(s);

            if (skipNaNInf && (!double.IsFinite(fs)))
                continue; // 或者 throw，看你调试策略

            double term = V[i] * fs;

            double y = term - c;
            double temp = sum + y;
            c = (temp - sum) - y;
            sum = temp;
        }

        return factor * sum;
    }

    // ===== 5) 向量反演（避免每次 Multiply 生成太多临时对象）=====
    public static Vector<double> InverseLaplaceVector(double t, int n, Func<double, Vector<double>> F, bool skipNaNInf = true)
    {
        if (t <= 0) throw new ArgumentException("t 必须 > 0");
        var V = GetVCoefficients(n);

        double ln2 = Math.Log(2.0);
        double factor = ln2 / t;

        Vector<double>? acc = null;

        for (int i = 1; i <= n; i++)
        {
            double s = i * factor;
            Vector<double> fi = F(s);

            if (acc == null)
            {
                acc = fi.Multiply(V[i]);
            }
            else
            {
                // 直接做 axpy：acc += V[i] * fi
                acc = acc.Add(fi.Multiply(V[i]));
            }
        }

        return acc!.Multiply(factor);
    }
}
