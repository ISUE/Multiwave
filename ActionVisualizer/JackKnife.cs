using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GestureTests.Data;
using GestureTests.Gesture;
using GestureTests.Util;
using GestureTests.Types;
using GestureTests.Experiment;
using System.Windows;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace ActionVisualizer
{
    class JackKnife
    {
        private List<Gesture> templates;
        private float pruned, tested;

        public JackKnife()
        {
            templates = new List<Gesture>();
            pruned = 0.0f;
            tested = 0.0f;
        }

        public void InitializeFromFolder(string path)
        {
            List<UserDataSet> alldata = DataLoader.LoadGestureDataFrom(path);
            foreach(UserDataSet ud in alldata)
            {
                foreach(GestureSample gs in ud.TrainingSamples)
                {
                    Add(gs.StrokePoints, gs.Gesture.ToString());
                }
                break;
            }
            Normalize();   
        }

        public void InitializeRawFromFolder(string path)
        {
            List<UserDataSet> alldata = DataLoader.LoadGestureDataFrom(path);
            foreach (UserDataSet ud in alldata)
            {
                foreach (GestureSample gs in ud.TrainingSamples)
                {
                    List<Vector<float>> data = new List<Vector<float>>();
                    foreach (var rd in gs.RawData)
                        data.Add(Vector<float>.Build.DenseOfArray(rd));

                    Add(data, gs.Gesture.ToString());
                }
                break;
            }
            Normalize();
        }
        
        public float CrossValidateDataset()
        {
            float errors = 0;
            for(int ii = 0; ii < templates.Count; ii++)
            {
                Gesture candidate = templates[ii];
                templates.RemoveAt(ii);

                var result = Classify(candidate);
                if (result.Item1.gname != candidate.gname)
                {
                    Console.WriteLine(result.Item1.gname + " " + candidate.gname);
                    errors++;
                }
                templates.Insert(ii, candidate);
            }
            return errors / (templates.Count);
        }   

        public void Add(Gesture alpha)
        {
            templates.Add(alpha);            
        }

        public void Add(List<Vector2> points, string label)
        {            
            Add(new Gesture(points, label));
        }

        public void Add(List<Point> temp, string label)
        {
            List<Vector2> points = new List<Vector2>();
            foreach (Point p in temp)
                points.Add(new Vector2((float) p.X, (float) p.Y));

            Add(new Gesture(points,label));
        }

        public void Add(List<Vector<float>> points, string label)
        {
            Add(new Gesture(points, label));
        }

        //Doesn't actually normalize but not going to refactor until Eugene does
        public void Normalize()
        {
            Random rand = new Random();
            if (templates.Count < 1)
                return;
            int n = templates[0].pts.Count() - 1;
            int m = (n + 1) / 2;

            for (int ii = 0; ii < 100; ii++)
            {
                List<Vector<float>> vecs = new List<Vector<float>>();
                n = templates[0].pts.Count() - 1;
                while (n > 0)
                {
                    int count = Math.Min(n, rand.Next(m, m));
                    //Get random element from templates
                    var t = templates.OrderBy(x => Guid.NewGuid()).FirstOrDefault();

                    var start = rand.Next(0, t.vecs.Count - count);
                    vecs.AddRange(t.vecs.GetRange(start, count));
                    n -= count;
                }
                foreach (var t in templates)
                {
                    float score = DTW_Distance(t.vecs, vecs, Gesture.r);
                    t.rejection.Add(score);
                }
            }

            foreach (var t in templates)
            {
                float mean = t.rejection.Average();
                float std = StdDev(t.rejection);
                t.rejection_threshold = mean - 0.0f * std;
            }
        }

        public Tuple<float, List<float>> UpperBound(Gesture candidate, Gesture template)
        {
            List<float> boundary = new List<float>();
            for (int ii = 0;
                ii < candidate.vecs.Count() &&
                ii < template.lower.Count() &&
                ii < template.upper.Count();
                ii++)
            {
                Vector<float> v1 = candidate.vecs[ii];
                Vector<float> lower = template.lower[ii];
                Vector<float> upper = template.upper[ii];

                var v2 = Vector<float>.Build.SameAs(v1);
                for (int jj = 0; jj < v2.Count(); jj++)
                {
                    if (v1[jj] > 0)
                        v2[jj] = upper[jj];
                    else
                        v2[jj] = lower[jj];
                }

                boundary.Add(1.0f - Math.Min(1.0f, v1.DotProduct(v2)));
            }
            return new Tuple<float, List<float>>(boundary.Sum(), boundary);
        }

        public Tuple<Gesture, float> Classify(Gesture candidate)
        {
            Gesture best = new Gesture(candidate.raw_pts, "");

            var temp_templates = templates;

            foreach (var t in temp_templates)
            {
                var out_tuple = UpperBound(candidate, t);
                t.bound = out_tuple.Item1;
                t.boundary = out_tuple.Item2;
            }

            temp_templates = temp_templates.OrderBy(x => x.bound).ToList();


            float best_score = float.PositiveInfinity;
            foreach (Gesture template in temp_templates)
            {
                tested += 1.0f;
                if (template.bound > best_score || template.bound > template.rejection_threshold)
                {
                    pruned += 1.0f;
                    continue;
                }

                float score = DTW_Distance(candidate.vecs, template.vecs, Gesture.r);

                if (score >= template.rejection_threshold)
                    continue;

                if (score < best_score)
                {
                    best = template;
                    best_score = score;
                }
            }
            return new Tuple<Gesture, float>(best, best_score);
        }

        public List<Tuple<Gesture, float>> Classify(Gesture candidate, int k)
        {

            var outList = new List<Tuple<Gesture, float>>();
            var temp_templates = templates;

            foreach (var t in temp_templates)
            {
                var out_tuple = UpperBound(candidate, t);
                t.bound = out_tuple.Item1;
                t.boundary = out_tuple.Item2;
            }

            temp_templates = temp_templates.OrderBy(x => x.bound).ToList();


            float best_score = float.PositiveInfinity;
            foreach (Gesture template in temp_templates)
            {
                tested += 1.0f;
                if (template.bound > best_score || template.bound > template.rejection_threshold)
                {
                    pruned += 1.0f;
                    continue;
                }

                float score = DTW_Distance(candidate.vecs, template.vecs, Gesture.r);

                if (score >= template.rejection_threshold)
                    continue;

                outList.Add(new Tuple<Gesture, float>(template, score));
            }
            outList.Sort((x, y) => x.Item2.CompareTo(y.Item2));
            return outList.GetRange(0, Math.Min(k,outList.Count()));
        }


        public float DTW_Distance(List<Vector<float>> candidate, List<Vector<float>> template, int r)
        {
            int m = candidate.Count + 1, n = candidate.Count + 1;
            float[,] dtw = new float[m, n];

            for (int ii = 0; ii < m; ii++)
                for (int jj = 0; jj < n; jj++)
                    dtw[ii, jj] = float.PositiveInfinity;
            dtw[0, 0] = 0;
            for (int ii = 1; ii < m; ii++)
            {
                for (int jj = Math.Max(1, ii - r); jj < Math.Min(m, ii + r); jj++)
                {
                    float cost = 1-candidate[ii - 1].DotProduct(template[jj - 1]);
                    float min = Math.Min(dtw[ii - 1, jj], Math.Min(dtw[ii, jj - 1], dtw[ii - 1, jj - 1]));
                    dtw[ii, jj] = cost + min;
                }
            }

            return dtw[candidate.Count, template.Count];
        }

        public float DTW_Distance_Max(List<Vector<float>> candidate, List<Vector<float>> template, int r)
        {
            int m = candidate.Count + 1, n = candidate.Count + 1;
            float[,] dtw = new float[m, n];

            for (int ii = 0; ii < m; ii++)
                for (int jj = 0; jj < n; jj++)
                    dtw[ii, jj] = float.NegativeInfinity;
            dtw[0, 0] = 0;
            for (int ii = 1; ii < m; ii++)
            {
                for (int jj = Math.Max(1, ii - r); jj < Math.Min(m, ii + r); jj++)
                {
                    float cost = candidate[ii - 1].DotProduct(template[jj - 1]);
                    float min = Math.Max(dtw[ii - 1, jj], Math.Max(dtw[ii, jj - 1], dtw[ii - 1, jj - 1]));
                    dtw[ii, jj] = cost + min;
                }
            }

            return dtw[candidate.Count, template.Count];
        }

        private float StdDev(IEnumerable<float> values)
        {
            float ret = 0;
            if (values.Count() > 0)
            {
                //Compute the Average      
                float avg = values.Average();
                //Perform the Sum of (value-avg)_2_2      
                float sum = (float)values.Sum(d => Math.Pow(d - avg, 2));
                //Put it all together      
                ret = (float)Math.Sqrt((sum) / (values.Count() - 1));
            }
            return ret;
        }
    }
    
    class Gesture
    {        
        /*
        public static int resample_cnt = 20;
        public static int r = 4;//(resample_cnt / 10);
        */
        
        public static int resample_cnt = 16;
        public static int r = (resample_cnt / 10);
        
        public string gname { get; set; }
        public List<Vector<float>> raw_pts;
        public List<Vector<float>> pts;
        public List<Vector<float>> vecs;
        public List<Vector<float>> lower, upper;

        public List<float> rejection;
        public float rejection_threshold;
        public float bound;
        internal List<float> boundary;

        public Gesture(List<Vector2> points, string label)
        {
            var temp = new List<Vector<float>>();
            foreach (Vector2 p in points)
            {
                float[] data = { p.X, p.Y };
                temp.Add(Vector<float>.Build.Dense(data));
            }
            gname = label;
            raw_pts = temp;
            pts = Resample(temp, resample_cnt);
            vecs = new List<Vector<float>>();
            for (int ii = 1; ii < pts.Count; ii++)
            {
                vecs.Add(pts[ii].Subtract(pts[ii - 1]).Normalize(2));
            }
            Tuple<List<Vector<float>>, List<Vector<float>>> out_tuple = Envelop(vecs, r);
            lower = out_tuple.Item1;
            upper = out_tuple.Item2;

            rejection = new List<float>();
        }

        public Gesture(List<Point> points, string label)
        {
            var temp = new List<Vector<float>>();
            foreach (Point p in points)
            {
                float[] data = { (float)p.X, (float)p.Y };
                temp.Add(Vector<float>.Build.Dense(data));
            }
            gname = label;
            raw_pts = temp;
            pts = Resample(temp, resample_cnt);
            vecs = new List<Vector<float>>();
            for (int ii = 1; ii < pts.Count; ii++)
            {
                vecs.Add(pts[ii].Subtract(pts[ii - 1]).Normalize(2));
            }
            Tuple<List<Vector<float>>, List<Vector<float>>> out_tuple = Envelop(vecs, r);
            lower = out_tuple.Item1;
            upper = out_tuple.Item2;

            rejection = new List<float>();
        }

        public Gesture(List<Vector<float>> temp, string label)
        {
            gname = label;
            raw_pts = temp;
            pts = Resample(temp, resample_cnt);
            vecs = new List<Vector<float>>();
            for (int ii = 1; ii < pts.Count; ii++)
            {
                vecs.Add(pts[ii].Subtract(pts[ii - 1]).Normalize(2));
            }
            Tuple<List<Vector<float>>, List<Vector<float>>> out_tuple = Envelop(vecs, r);
            lower = out_tuple.Item1;
            upper = out_tuple.Item2;

            rejection = new List<float>();
        }

        private float PathLength(List<Vector<float>> points)
        {
            float length = 0;
            for (int ii = 1; ii < points.Count; ii++)
            {
                length += (float)points[ii].Subtract(points[ii - 1]).L2Norm();
            }
            return length;
        }

        private List<Vector<float>> Resample(List<Vector<float>> temp, int n)
        {
            List<Vector<float>> points = new List<Vector<float>>();
            foreach (Vector<float> v in temp)
                points.Add(Vector<float>.Build.DenseOfVector(v));
            List<Vector<float>> ret = new List<Vector<float>>();
            ret.Add(Vector<float>.Build.DenseOfVector(points[0]));
            float I = PathLength(points) / (n - 1.0f);          
            float D = 0.0f;
            int ii = 1;
            while (ii < points.Count && I > 0)
            {
                float d = (float)points[ii].Subtract(points[ii - 1]).L2Norm();
                if (D + d >= I)
                {
                    Vector<float> vec = points[ii].Subtract(points[ii - 1]);
                    float t = (I - D) / d;

                    if (float.IsNaN(t))
                        t = 0.5f;

                    Vector<float> q = points[ii - 1] + t * vec;
                    ret.Add(q);
                    points.Insert(ii, q);
                    D = 0.0f;
                }
                else
                {
                    D += d;
                }
                ++ii;
            }
            while (ret.Count < n)
                ret.Add(points.Last());
            return ret;
        }

        private Tuple<List<Vector<float>>, List<Vector<float>>> Envelop(List<Vector<float>> vecs, int r)
        {
            int n = vecs.Count;
            int m = vecs[0].Count;

            List<Vector<float>> upper = new List<Vector<float>>();
            List<Vector<float>> lower = new List<Vector<float>>();

            for (int ii = 0; ii < n; ii++)
            {
                Vector<float> maximum = Vector<float>.Build.Dense(m, float.NegativeInfinity);
                Vector<float> minimum = Vector<float>.Build.Dense(m, float.PositiveInfinity);
                for (int jj = Math.Max(0, ii - r); jj < Math.Min(ii + r + 1, n); jj++)
                {
                    maximum = Maximum(maximum, vecs[jj]);
                    minimum = Minimum(minimum, vecs[jj]);
                }
                upper.Add(maximum);
                lower.Add(minimum);
            }

            return new Tuple<List<Vector<float>>, List<Vector<float>>>(lower, upper);
        }

        private Vector<float> Maximum(Vector<float> a, Vector<float> b)
        {
            Vector<float> max = Vector<float>.Build.SameAs(a);
            for (int ii = 0; ii < a.Count(); ii++)
                max[ii] = Math.Max(a[ii], b[ii]);
            return max;
        }
        private Vector<float> Minimum(Vector<float> a, Vector<float> b)
        {
            Vector<float> min = Vector<float>.Build.SameAs(a);
            for (int ii = 0; ii < a.Count(); ii++)
                min[ii] = Math.Min(a[ii], b[ii]);
            return min;
        }

        public int CompareTo(Gesture other)
        {
            return bound.CompareTo(other.bound);
        }
    }
}
