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
            ReservoirConfigDoubleMedia[] configs = new ReservoirConfigDoubleMedia[2];
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

            ReservoirConfigDoubleMedia config1 = new ReservoirConfigDoubleMedia(h: 6, k: 0.001, Bo: 1.2, mu: 12, ctm: 0.6e-4, ctf: 4e-4, Cs: 0.01, rw: 0.12, re: 20000, phiM: 0.1, phiF: 0.014, s: 0.5, sf: 0.5, xf: 100)
            {
                Pi = 32
            };
            ReservoirConfigDoubleMedia config2 = new ReservoirConfigDoubleMedia(h: 6, k: 0.0001, 1.2, 15, ctm: 0.5e-4, ctf: 3e-4, Cs: 0.01, 0.12, 20000, phiM: 0.06, phiF: 0.01, s: 0.5, sf: 0.5, xf: 100)
            {
                Pi = 32
            }
            ;//4, 0.0002, 2.5618E-4, 12, 5e-4, 1.4, 0.1, 0.12, 1.5
            
            config1.factorlessness(config2, omega: 0.7, lambdaF: 1e-2,zeta:0.7);
            config2.factorlessness(config1, omega: 0.3, lambdaF: 1e-2, zeta:0.3);
            configs[0] = config1;
            configs[1] = config2;
            //double[,] qdbint = new double[2, 161];
            for (int i = -40; i <= 160; i++)
            {
                double tD = Math.Pow(10, 0.1 * i) * (3.6 * config1.eta) / config1.xf / config1.xf;
                //Debug.WriteLine("t:" + Math.Pow(10, 0.1 * i) + "   无因次化：" + (3.6 * config1.kf /config1.mu/(config1.totalPhiCt)/ config1.xf / config1.xf));
                double Pwf = MathUnit.Getf(Math.Pow(10, 0.1 * i) * (3.6 * config1.eta) / config1.xf / config1.xf, 6, configs);
                double DPwf = MathUnit.GetDf(Pwf, Math.Pow(10, 0.1 * i) * (3.6 * config1.eta / config1.xf / config1.xf), 6, configs);
                
                
                double[] qwf = new double[2];
                qwf = MathUnit.GetCQ(Math.Pow(10, 0.1 * i) * (3.6 * config1.eta) / config1.xf / config1.xf, 6, configs);
                for (int j = 0; j < 2; j++)
                {
                    qwf[j] = qwf[j] * (-1);
                }

                double tDd = 2*tD/(config1.reD* config1.reD - 1)/(Math.Log(config1.reD)-0.5);
                var q = MathUnit.产量递减GetQLapalce(tDd, configs);
                /*qdbint[0, i + 40] = q[0];
                qdbint[1, i + 40] = q[1];*/
                var qdi = MathUnit.产量递减GetQDi(q, tDd);
                var dq = MathUnit.产能递减GetdQDi(q, tDd, configs);
                //Console.WriteLine(i);
                //textBox7.Text += Math.Pow(10, 0.1 * i) + "\t" + Pwf + "\t" + DPwf + "\t" + Qwf + "\t" + DQwf + "\r\n";//
                textBox1.Text += tD + "\t" + Pwf + "\t" + DPwf + "\t" + tDd + "\t" + q[0] + "\t" + q[1] + "\t" + dq[0] + "\t" + dq[1] + "\r\n";//

            }
        }
    }
}
