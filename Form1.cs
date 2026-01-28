using ShaleOilWellTest;

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
            ReservoirConfigDoubleMedia[] configs = InitLayers();
            foreach (var config in configs)
            {
                config.factorlessness(configs);
            }
            var model = SelectModel(checkBox1.Checked, configs);
            var sb = new System.Text.StringBuilder();


            //double[,] qdbint = new double[2, 161];
            for (int i = -80; i <= 80; i++)
            {
                double tD = Math.Pow(10, 0.1 * i) * (3.6 * configs[0].eta) / configs[0].rw / configs[0].rw;
                AppendResult(sb, tD, model.pwf, model.derivative);
            }
            textBox1.Text = sb.ToString();
        }

        private static ReservoirConfigDoubleMedia[] InitLayers()
        {
            var list = new List<ReservoirConfigDoubleMedia>();
            double theta = double.Pi / 6;


            list.Add(new ReservoirConfigDoubleMedia(index:0,h: 8, k: 0.28, Bo: 1.2, mu: 3.5, ctm: 8e-4, ctf: 2e-3, Cs: 0.0001, rw: 0.1, re: 400, phiM: 0.1104, phiF: 0.032, s: 5, sf: 1, xf: 100)
            {
                theta = theta,
                lambdaF = 1e-4,
                Pi = 32
            });
            list.Add(new ReservoirConfigDoubleMedia(index: 1,h: 8, k: 0.60, Bo: 1.2, mu: 3.5, ctm: 6e-4, ctf: 1e-3, Cs: 0.0001, rw: 0.1, re: 400, phiM: 0.09, phiF: 0.02, s: 5, sf: 0.2, xf: 30)
            {
                theta = theta,
                lambdaF = 1e-4,
                Pi = 36
            });
            list.Add(new ReservoirConfigDoubleMedia(index: 2, h: 5, k: 0.04, Bo: 1.2, mu: 3.5, ctm: 2e-4, ctf: 0.7e-3, Cs: 0.0001, rw: 0.1, re: 400, phiM: 0.03, phiF: 0.004, s: 5, sf: 0.1, xf: 30)
            {
                theta = theta,
                lambdaF = 1e-4,
                Pi = 40
            });
            return list.ToArray();
        }

        private static void AppendResult(System.Text.StringBuilder sb, double tD, Func<double, double> pwfFunc, Func<double, double, double> derivativeFunc)
        {
            double pwfD = pwfFunc(tD);
            double dpwfD = derivativeFunc(pwfD, tD);
            sb.Append($"{tD}\t{pwfD}\t{dpwfD}\r\n");
        }

        /// <summary>
        /// 选择计算模型。后续添加新模型时只需在此扩展，不必修改主循环。
        /// </summary>
        private static (Func<double, double> pwf, Func<double, double, double> derivative) SelectModel(bool isMultiLayer, ReservoirConfigDoubleMedia[] configs)
        {
            const int stehfestN = 10;
            if (isMultiLayer)
            {
                double Pwf(double t) => MathUnit.MultiTestGetf(t, stehfestN, configs);
                double Dpwf(double pwf, double t) => MathUnit.MultiTestGetDf(pwf, t, stehfestN, configs);
                return (Pwf, Dpwf);
            }
            else
            {
                double Pwf(double t) => MathUnit.Getf_cui(t, stehfestN, configs[0]);
                double Dpwf(double pwf, double t) => MathUnit.GetDf_cui(pwf, t, stehfestN, configs[0]);
                return (Pwf, Dpwf);
            }
        }
    }
}
