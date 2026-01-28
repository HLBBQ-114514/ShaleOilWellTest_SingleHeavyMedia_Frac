using MathNet.Numerics;
using ShaleOilWellTest;
using System.Diagnostics;

namespace ShaleOilWellTest
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            checkBox1.Text = "是否多层";

        }

        private void button1_Click(object sender, EventArgs e)
        {
            textBox1.Clear();
            ReservoirConfigDoubleMedia[] configs = null!;
            configs = init();
            foreach (var config in configs)
            {
                config.factorlessness(configs);
            }
            //小层初始化
            static ReservoirConfigDoubleMedia[] init() 
            {
                var list = new List<ReservoirConfigDoubleMedia>();
                double theta = double.Pi/6;
                
                ReservoirConfigDoubleMedia.Pi = 32;
                //(h: 2, k: 0.06, Bo: 1.2, mu: 3.5, ctm: 0.5e-4, ctf: 3e-3, Cs: 0.01, rw: 0.12, re: 3000, phiM: 0.13, phiF: 0.002, s: 0, sf: 0.1, xf: 40)
                list.Add(new ReservoirConfigDoubleMedia(h: 8, k: 0.28, Bo: 1.2, mu: 3.5, ctm: 8e-4, ctf: 2e-3, Cs: 1e-6, rw: 0.1, re: 800, phiM: 0.1104, phiF: 0.032, s: 0.1, sf: 1, xf: 100)
                {
                    /*                omega = omega1,
                                    zeta = 0.1,*/
                    theta = theta,
                    lambdaF = 1e-4,
                });
                //(h: 1.2, k: 0.4, Bo: 1.2, mu: 3.5, ctm: 0.9e-4, ctf: 4e-3, Cs: 0.01, rw: 0.12, re: 3000, phiM: 0.1304, phiF: 0.062, s: -0.3, sf: 0.1, xf: 100)
                list.Add(new ReservoirConfigDoubleMedia(h: 8, k: 0.60, Bo: 1.2, mu: 3.5, ctm: 6e-4, ctf: 1e-3, Cs: 0.0001, rw: 0.1, re: 800, phiM: 0.09, phiF: 0.02, s: 1, sf: 0.2, xf: 30)
                {
                    /*                omega = omega2,
                                    zeta = 0.9,*/
                    theta = theta,
                    lambdaF = 1e-4,

                });
                list.Add(new ReservoirConfigDoubleMedia(h: 5, k: 0.04, Bo: 1.2, mu: 3.5, ctm: 2e-4, ctf: 0.7e-3, Cs: 0.0001, rw: 0.1, re: 800, phiM: 0.03, phiF: 0.004, s: 1, sf: 0.1, xf: 30)
                {
                    theta = theta,
                    lambdaF = 1e-4
                });
                //4, 0.0002, 2.5618E-4, 12, 5e-4, 1.4, 0.1, 0.12, 1.5
                return list.ToArray();
            }
            double omega1 = 0.5;
            double omega2 = (1 - omega1);
            

            // === 新增产量递减变量 ===
            List<double> tList = new List<double>();
            List<double> qDdList = new List<double>();
            List<double> qDdiList = new List<double>();
            List<double> qDdidList = new List<double>();
            List<double> qDcumList = new List<double>() { 0 }; // 累计无因次产量初始值为0
            //double integral = 0; // 累计积分项
            var sb = new System.Text.StringBuilder();


            //double[,] qdbint = new double[2, 161];
            for (int i = -100; i <= 40; i++)
            {
                double tD = Math.Pow(10, 0.1 * i) * (3.6 * configs[0].eta) / configs[0].rw / configs[0].rw;
                #region 多层合采试井计算
                /*
                //Debug.WriteLine("t:" + Math.Pow(10, 0.1 * i) + "   无因次化：" + (3.6 * configs[0].kf /configs[0].mu/(configs[0].totalPhiCt)/ configs[0].xf / configs[0].xf));
                double PwfD = MathUnit.Getf(tD, 6, configs);                
                if (double.IsNaN(PwfD)) continue;
                double DPwfD = MathUnit.GetDf(PwfD, tD, 6, configs);

                
                double[] qwf = new double[2];
                qwf = MathUnit.GetCQ(tD, 6, configs);
                for (int j = 0; j < 2; j++)
                {
                    qwf[j] = qwf[j] * (-1);
                }

                //double tDd = 2*tD/(configs[0].reD* configs[0].reD - 1)/(Math.Log(configs[0].reD)-0.5);
                
                var qD = MathUnit.产量递减GetQLapalce(tD, configs);
                if(qD.Any(x => double.IsNaN(x))) continue;
                double[] q = new double[qD.Length];
                double qDd = 0;
                for (int j = 0; j < qD.Length; j++)
                {
                    q[j] = qD[j] * 2 * double.Pi * configs[0].avgLambda * (ReservoirConfigDoubleMedia.Pi-5);
                    qDd += q[j];
                }
                qDdList.Add(qDd);
                tList.Add(tD);

                // 2️⃣ 累计无因次产量（积分平均）
                if (qDdList.Count == 1)
                {
                    qDdiList.Add(0);
                }
                else
                {
                    int k = qDdList.Count - 1;
                    double dt = tList[k] - tList[k - 1];
                    integral += 0.5 * (qDdList[k] + qDdList[k - 1]) * dt;
                    double qDdi = integral / tList[k];
                    qDdiList.Add(qDdi);
                    double qDcum = qDcumList[^1] + 0.5 * (qDdList[k] + qDdList[k - 1]) * dt;
                    qDcumList.Add(qDcum);
                }

                // 3️⃣ 累计导数（Blasingame 理论关系）
                // q_Ddid = q_Dd - q_Ddi
                double qDdid = qDdiList[^1] - qDd;
                qDdidList.Add(qDdid);
                var s = $"{qDd}\t{qDdiList[^1]}\t{qDdid}";//\t{qDcumList[^1]}

                var qj = MathUnit.GetQ(tD, 6, configs);
                qj.ToString();
                //var sq = string.
                var sq = string.Join("\t ", qj);
                *//*qdbint[0, i + 40] = q[0];
                qdbint[1, i + 40] = q[1];*//*
                //var qdi = MathUnit.产量递减GetQDi(qD, tD);
                //var dq = MathUnit.产能递减GetdQDi(qD, tD, configs);
                //Console.WriteLine(i);
                //textBox7.Text += Math.Pow(10, 0.1 * i) + "\t" + Pwf + "\t" + DPwf + "\t" + Qwf + "\t" + DQwf + "\r\n";//
                //使用lambda语句法输出qD 并且使用math.round截取小数点
                //string qtext = $"{Math.Round(q[0], 4)} \t  {Math.Round(q[1], 4)} \t {Math.Round(q[2], 4)}";//
                textBox1.Text += tD + "\t" + Math.Abs(Math.Round(PwfD,12)) + "\t" + Math.Abs(Math.Round(DPwfD,12)) + "\t" + s + "\t" + sq+"\r\n";//
               */
                #endregion
                #region cui调试代码
                if(checkBox1.Checked == false)
                {
                    double pwfD_debug = MathUnit.Getf_cui(tD, 10, configs[0]);
                    double DpwfD_debug = MathUnit.GetDf_cui(pwfD_debug, tD, 10, configs[0]);
                    sb.Append($"{tD}\t{pwfD_debug}\t{DpwfD_debug}\n");

                }
                else
                {
                    double pwfD_multi_debug = MathUnit.MultiTestGetf(tD, 10, configs);
                    double DpwfD_multi_debug = MathUnit.MultiTestGetDf(pwfD_multi_debug, tD, 10, configs);
                    sb.Append($"{tD}\t{pwfD_multi_debug}\t{DpwfD_multi_debug}\r\n");
                }

                #endregion
            }
            textBox1.Text = sb.ToString();
        }
    }
}
