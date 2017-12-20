using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using InoueLab;
using System.IO;
using System.Windows.Forms.DataVisualization.Charting;

namespace WindowsFormsApplication9
    {
    public partial class Form1 : Form
        {
        public Form1()
            {




            InitializeComponent();
            double[] _X = new double[CONST.N];
            double[,] _designmatrix = new double[CONST.N, CONST.M];
            double[] _t = new double[CONST.N];
            double[,] _weight = new double[CONST.M, 1];
            MAKEDATA(ref _designmatrix, ref _t, ref _X);
            _weight = Calc(_designmatrix, _t);
            Plot(_X, _t, _weight);
            }
        private void Plot(double[] _X, double[] _t, double[,] w)
            {
            chart1.Series.Clear();
            chart1.Series.Add("sin");
            chart1.Series.Add("t");
            chart1.Series.Add("least square method");

            chart1.Series["sin"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            chart1.Series["t"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Point;
            chart1.Series["least square method"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;


            for (double i = 0; i <= 2 * Math.PI; i += Math.PI / 360)
                {
                chart1.Series["sin"].Points.AddXY(i, Math.Sin(i));
                double x = i;
                double y = w[0, 0];
                for (int j = 1; j < CONST.M; j++)
                    {
                    y += x * w[j, 0];
                    x *= i;
                    }
                chart1.Series["least square method"].Points.AddXY(i, y);

                }
            for (int i = 0; i < CONST.N; i++)
                {
                chart1.Series["t"].Points.AddXY(_X[i], _t[i]);

                }
            }
        private static class CONST
            {
            public const int N = 10;
            public const int M = 4;
            }

        private class RandomBoxMuller
            {
            private Random random;

            public RandomBoxMuller()
                {
                random = new Random(Environment.TickCount);
                }

            public RandomBoxMuller(int seed)
                {
                random = new Random(seed);
                }

            public double next(double mu = 0.0, double sigma = 1.0, bool getCos = true)
                {
                if (getCos)
                    {
                    double rand = 0.0;
                    while ((rand = random.NextDouble()) == 0.0) ;
                    double rand2 = random.NextDouble();
                    double normrand = Math.Sqrt(-2.0 * Math.Log(rand)) * Math.Cos(2.0 * Math.PI * rand2);
                    normrand = normrand * sigma + mu;
                    return normrand;
                    }
                else
                    {
                    double rand;
                    while ((rand = random.NextDouble()) == 0.0) ;
                    double rand2 = random.NextDouble();
                    double normrand = Math.Sqrt(-2.0 * Math.Log(rand)) * Math.Sin(2.0 * Math.PI * rand2);
                    normrand = normrand * sigma + mu;
                    return normrand;
                    }
                }

            public double[] nextPair(double mu = 0.0, double sigma = 1.0)
                {
                double[] normrand = new double[2];
                double rand = 0.0;
                while ((rand = random.NextDouble()) == 0.0) ;
                double rand2 = random.NextDouble();
                normrand[0] = Math.Sqrt(-2.0 * Math.Log(rand)) * Math.Cos(2.0 * Math.PI * rand2);
                normrand[0] = normrand[0] * sigma + mu;
                normrand[1] = Math.Sqrt(-2.0 * Math.Log(rand)) * Math.Sin(2.0 * Math.PI * rand2);
                normrand[1] = normrand[1] * sigma + mu;
                return normrand;
                }

            }
        private void MAKEDATA(ref double[,] _designmatrix, ref double[] _t, ref double[] _X)
            {
            var r = new RandomBoxMuller();
            RandomMT rand = new RandomMT(1);
            for (int i = 0; i < CONST.N; i++)
                {
                _X[i] = rand.Double() * Math.PI * 2;
                _t[i] = System.Math.Sin(_X[i]) + r.next();
                }
            for (int j = 0; j < CONST.N; j++)
                {
                for (int i = 0; i < CONST.M; i++)
                    {
                    _designmatrix[j, i] = System.Math.Pow(_X[j], i);
                    if (i == 0)
                        {
                        _designmatrix[j, i] = 1;
                        }
                    }
                }

            }

        private double[,] Calc(double[,] _designmatrix, double[] _t)
            {
            double[,] _tT = new double[CONST.N, 1];
            double[] _y = new double[CONST.N];
            double[] _design = new double[CONST.M];
            double[,] _designmatrixT = new double[CONST.M, CONST.N];
            double[,] _I = new double[CONST.M, CONST.M];

            for (int j = 0; j < CONST.N; j++)
                {
                for (int i = 0; i < CONST.M; i++)
                    {
                    _designmatrixT[i, j] = _designmatrix[j, i];
                    }
                }

            for (int i = 0; i < CONST.N; i++)
                {
                _tT[i, 0] = _t[i];
                }

            double[,] _Designsquare = new double[CONST.M, CONST.M];
            double integer;

            for (int i = 0; i < CONST.M; i++)
                {
                for (int j = 0; j < CONST.M; j++)
                    {
                        integer = 0;
                        for (int n = 0; n < CONST.N; n++)
                            {
                            integer += _designmatrixT[i, n] * _designmatrix[n, j];
                            }

                    _Designsquare[i, j] = integer;
                    }

                }

            double buf;

            for (int i = 0; i < CONST.M; i++)
                {
                for (int j = 0; j < CONST.M; j++)
                    {
                    _I[i, j] = 0;
                    if (i == j) _I[i, j] = 1;
                    }
                }
            for (int i = 0; i < CONST.M; i++)
                {
                buf = 1 / _Designsquare[i, i];
                for (int j = 0; j < CONST.M; j++)
                    {
                    _Designsquare[i, j] *= buf;
                    _I[i, j] *= buf;
                    }
                for (int j = 0; j < CONST.M; j++)
                    {
                    if (i != j)
                        {
                        buf = _Designsquare[j, i];
                        for (int k = 0; k < CONST.M; k++)
                            {
                            _Designsquare[j, k] -= _Designsquare[i, k] * buf;
                            _I[j, k] -= _I[i, k] * buf;
                            }
                        }
                    }
                }

            double[,] _design_t = new double[CONST.M, 1];
            double _integer2;
            for (int i = 0; i < CONST.M; i++)
                {
                _integer2 = 0;

                for (int j = 0; j < CONST.N; j++)
                    {
                    _integer2 += _designmatrixT[i, j] * _tT[j, 0];
                    }
                _design_t[i, 0] = _integer2;
                }


            double[,] _weight = new double[CONST.M, 1];

            for (int i = 0; i < CONST.M; i++)
                {
                _integer2 = 0;
                for (int j = 0; j < CONST.M; j++)
                    {
                    _integer2 += _I[i, j] * _design_t[j, 0];
                    }
                _weight[i, 0] = _integer2;
                }

            return _weight;
            }

        }
    }



