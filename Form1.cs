using System.Diagnostics;

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
            
            ReservoirConfig[] configs = new ReservoirConfig[2];
            ReservoirConfig config1 = new(2, 0.0001, 1.5618E-4, 15, 3e-4, 0.12, 200, 0.1, 0.06, 2, 0.1, 100)
            {
                Pi = 32
            };
            ReservoirConfig config2 = new(4, 0.0002, 2.5618E-4, 12, 5e-4, 0.12, 200, 0.1, 0.1, 1.5, 0.1, 100)
            {
                Pi = 32
            };//4, 0.0002, 2.5618E-4, 12, 5e-4, 1.4, 0.1, 0.12, 1.5
            config1.factorlessness(config2);
            config2.factorlessness(config1);
            configs[0] = config1;
            configs[1] = config2;

            for (int i = -40; i <= 70; i++)
            {
                Debug.WriteLine("t:" + Math.Pow(10, 0.1 * i) + "   无因次化：" + (3.6 * config1.eta / config1.xf / config1.xf));
                double Pwf = MathUnit.Getf(Math.Pow(10, 0.1 * i) * (3.6 * config1.eta / config1.xf / config1.xf), 6, configs);
                double DPwf = MathUnit.GetDf(Pwf, Math.Pow(10, 0.1 * i) * (3.6 * config1.eta / config1.xf / config1.xf), 6, configs);
                
                
                double[] qwf = new double[2];
                qwf = MathUnit.GetCQ(Math.Pow(10, 0.1 * i) * (3.6 * config1.eta / config1.xf / config1.xf), 6, configs);
                for (int j = 0; j < 2; j++)
                {
                    qwf[j] = qwf[j] * (-1);
                }
                //Console.WriteLine(i);
                //textBox7.Text += Math.Pow(10, 0.1 * i) + "\t" + Pwf + "\t" + DPwf + "\t" + Qwf + "\t" + DQwf + "\r\n";//
                textBox1.Text += Math.Pow(10, 0.1 * i) + "\t" + Pwf + "\t" + DPwf+"\t" + "\r\n";//+ qwf[0] + "\t" + qwf[1]

            }
        }
    }
}
