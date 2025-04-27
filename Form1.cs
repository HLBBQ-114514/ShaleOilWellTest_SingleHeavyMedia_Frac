using System.Diagnostics;
using ShaleOilWellTest;

namespace ShaleOilWellTest
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            ReservoirConfigDoubleMedia[] configs = new ReservoirConfigDoubleMedia[3];
            textBox1.Text += "t\tpwf\tdpwf\tt\tq1\tq2\tq3\tdq1\tdq2\tdq3\r\n";

            //小层初始化
            static void init() 
            {
                int n = 2;
                
                for (int i = 0; i < n; i++)
                {   
                    //待完成
                    //configs[i] = new ReservoirConfig();
                }
            }

            configs[0] = new ReservoirConfigDoubleMedia(h: 2.5, k: 0.01, Bo: 1.2, mu: 15, ctm: 0.5e-4, ctf: 3e-4, Cs: 0.01, rw: 0.12, re: 2000, phiM: 0.03, phiF: 0.002, s: 5.9, sf: 5.9, xf: 100)
            {
                Pi = 32,
                omega = 0.1,
                zeta = 0.1,
                lambdaF = 1e-2,
                
            };
            configs[1] = new ReservoirConfigDoubleMedia(h: 1.5, k: 0.19, Bo: 1.2, mu: 12, ctm: 0.6e-4, ctf: 4e-4, Cs: 0.01, rw: 0.12, re: 2000, phiM: 0.1304, phiF: 0.062, s: 5.9, sf: 5.9, xf: 100)
            {
                Pi = 32,
                omega = 0.6,
                zeta = 0.6,
                lambdaF = 1e-2,
                
            };
            configs[2] = new ReservoirConfigDoubleMedia(h: 3, k: 0.01, Bo: 1.2, mu: 15, ctm: 0.5e-4, ctf: 3e-4, Cs: 0.01, rw: 0.12, re: 2000, phiM: 0.03, phiF: 0.002, s: 5.9, sf: 5.9, xf: 100)
            {
                Pi = 32,
                omega = 0.3,
                zeta = 0.3,
                lambdaF = 1e-2

            };
            ;//4, 0.0002, 2.5618E-4, 12, 5e-4, 1.4, 0.1, 0.12, 1.5
            foreach(var config in configs) {
                config.factorlessness(configs);
            }

            //double[,] qdbint = new double[2, 161];
            for (int i = -60; i <= 100; i++)
            {
                double tD = Math.Pow(10, 0.1 * i) * (3.6 * configs[0].eta) / configs[0].xf / configs[0].xf;
                //Debug.WriteLine("t:" + Math.Pow(10, 0.1 * i) + "   无因次化：" + (3.6 * configs[0].kf /configs[0].mu/(configs[0].totalPhiCt)/ configs[0].xf / configs[0].xf));
                double Pwf = MathUnit.Getf(Math.Pow(10, 0.1 * i) * (3.6 * configs[0].eta) / configs[0].xf / configs[0].xf, 6, configs);
                double DPwf = MathUnit.GetDf(Pwf, Math.Pow(10, 0.1 * i) * (3.6 * configs[0].eta / configs[0].xf / configs[0].xf), 6, configs);
                
                
                double[] qwf = new double[2];
                qwf = MathUnit.GetCQ(Math.Pow(10, 0.1 * i) * (3.6 * configs[0].eta) / configs[0].xf / configs[0].xf, 6, configs);
                for (int j = 0; j < 2; j++)
                {
                    qwf[j] = qwf[j] * (-1);
                }

                //double tDd = 2*tD/(configs[0].reD* configs[0].reD - 1)/(Math.Log(configs[0].reD)-0.5);
                var q = MathUnit.产量递减GetQLapalce(tD, configs);
                /*qdbint[0, i + 40] = q[0];
                qdbint[1, i + 40] = q[1];*/
                var qdi = MathUnit.产量递减GetQDi(q, tD);
                var dq = MathUnit.产能递减GetdQDi(q, tD, configs);
                //Console.WriteLine(i);
                //textBox7.Text += Math.Pow(10, 0.1 * i) + "\t" + Pwf + "\t" + DPwf + "\t" + Qwf + "\t" + DQwf + "\r\n";//
                textBox1.Text += tD + "\t" + Pwf + "\t" + DPwf + "\t" + tD + "\t" + q[0] + "\t" + q[1] + "\t" + dq[0] + "\t" + dq[1] + "\r\n";//

            }
        }
    }
}
