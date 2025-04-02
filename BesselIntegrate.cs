using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.Integration;
using ShaleOilWellTest;

namespace ShaleOilWellTest_单重介质压裂封闭边界
{
    internal class BesselIntegrate
    {   
    /*
     * <summary>
     * 
     * </summary>
     * <param name="func">被积函数</param>
     * <param name="xDjStart">积分下限</param>
     * <param name="xDjEnd">积分上限</param>
     */
        public static double BesselInt(Func<double,double> func, double xDjStart,double xDjEnd)
        {
            return GaussLegendreRule.Integrate(func, xDjStart, xDjEnd, 32);
        }

    }
}
