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
using System.Threading;

namespace ActionVisualizer
{

    class RecognitionResult
    {
        public float tests;
        public float pruned_by_best_score;
        public float pruned_by_rejection;
        public Gesture template;
        public float score;

        public RecognitionResult()
        {
            tests = 0.0f;
            pruned_by_best_score = 0.0f;
            pruned_by_rejection = 0.0f;
            template = null;
            score = float.PositiveInfinity;
        }
    }


    class JackKnife
    {
        public List<Gesture> templates;

        public JackKnife()
        {
            templates = new List<Gesture>();           
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
            }
            Normalize();
        }
        public void InitializeFromSingleUser(string path, string uname)
        {
            List<UserDataSet> alldata = DataLoader.LoadGestureDataFrom(path);
            foreach (UserDataSet ud in alldata)
            {
                if (uname == "" || ud.Path.Contains(uname))
                {
                    foreach (GestureSample gs in ud.TrainingSamples)
                    {

                        List<Vector<float>> data = new List<Vector<float>>();
                        foreach (var rd in gs.RawData)
                            data.Add(Vector<float>.Build.DenseOfArray(rd));

                        Add(data, gs.Gesture.ToString());
                    }
                }
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
                if(result.template == null)
                {
                    Console.WriteLine("Nothing " + candidate.gname);
                    errors++;
                }
                else if (result.template.gname != candidate.gname)
                {
                    Console.WriteLine(result.template.gname + " " + candidate.gname);
                    errors++;
                }
                templates.Insert(ii, candidate);
            }
            return errors / (templates.Count);
        }

        public float CompareDatasets()
        {
            float errors = 0, errors_2D = 0, errors_3D = 0;
            List<Gesture> original = templates;
            List<Gesture> only2D = original.FindAll(x => x.Is2D() == true);
            List<Gesture> only3D = original.FindAll(x => x.Is2D() == false);

            Console.WriteLine();
            Console.WriteLine("Comparing 2D");

            templates = only2D;
            for (int ii = 0; ii < templates.Count; ii++)
            {
                Gesture candidate = templates[ii];
                templates.RemoveAt(ii);

                var result = Classify(candidate);
                if (result.template == null)
                {
                    Console.WriteLine("Nothing " + candidate.gname);
                    errors_2D++;
                }
                else if (result.template.gname != candidate.gname)
                {
                    Console.WriteLine(result.template.gname + " " + candidate.gname);
                    errors_2D++;
                }
                templates.Insert(ii, candidate);
            }

            Console.WriteLine();
            Console.WriteLine("Comparing 3D");

            templates = only3D;
            for (int ii = 0; ii < templates.Count; ii++)
            {
                Gesture candidate = templates[ii];
                templates.RemoveAt(ii);

                var result = Classify(candidate);
                if (result.template == null)
                {
                    Console.WriteLine("Nothing " + candidate.gname);
                    errors_3D++;
                }
                else if (result.template.gname != candidate.gname)
                {
                    Console.WriteLine(result.template.gname + " " + candidate.gname);
                    errors_3D++;
                }
                templates.Insert(ii, candidate);
            }

            templates = original;

            Console.WriteLine();
            Console.WriteLine("Comparing All");

            for (int ii = 0; ii < templates.Count; ii++)
            {
                Gesture candidate = templates[ii];
                templates.RemoveAt(ii);

                var result = Classify(candidate);
                if (result.template == null)
                {
                    Console.WriteLine("Nothing " + candidate.gname);
                    errors++;
                }
                else if (result.template.gname != candidate.gname)
                {
                    Console.WriteLine(result.template.gname + " " + candidate.gname);
                    errors++;
                }
                templates.Insert(ii, candidate);
            }

            Console.WriteLine();
            Console.WriteLine("Error Rates:");

            errors = errors / templates.Count;
            errors_2D = errors_2D / only2D.Count;
            errors_3D = errors_3D / only3D.Count;

            Console.WriteLine("2D: \t" + errors_2D);
            Console.WriteLine("3D: \t" + errors_3D);
            Console.WriteLine("All:\t" + errors);

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
                t.rejection_threshold = mean - 1.5f * std;
            }
        }

        public float LowerBound(
             Gesture candidate,
             Gesture template)
        {
            int m = candidate.vecs[0].Count;
            float ret = 0.0f;

            for (int ii = 0;
                 ii < candidate.vecs.Count();
                 ii++)
            {
                float tmp = 0.0f;

                for (int jj = 0;
                     jj < m;
                     jj++)
                {
                    if (candidate.vecs[ii][jj] >= 0.0f)
                        tmp += candidate.vecs[ii][jj] * template.upper[ii][jj];
                    else
                        tmp += candidate.vecs[ii][jj] * template.lower[ii][jj];
                }

                tmp = Clamp(tmp, -1.0f, 1.0f);
                ret += 1.0f - tmp;
            }

            return ret;
        }


        public RecognitionResult Classify(Gesture candidate)
        {
            int feature_cnt = candidate.features.Count();
            //
            // calculate lower bounds for each template based on envelop
            // and sort the list based on the boundaries  
            //
            var temp_templates = templates;
            foreach (var t in temp_templates)
            {
                float score = 1.0f;
                for (int ii = 0;
                    ii < feature_cnt;
                    ii++)
                {
                    Vector<float> fvec = candidate.features[ii];
                    score *= 1.0f / fvec.DotProduct(t.features[ii]);
                }
                t.fscore = score;
                t.lower_bound = LowerBound(candidate, t);
                //Debug.Log(t.gname + ":" + t.lower_bound.ToString() + ", " + t.rejection_threshold.ToString() ); 
            }

            //Debug.Log("\n");

            temp_templates = temp_templates.OrderBy(x => x.lower_bound).ToList();

            //
            // now fully evaluate each template
            //
            RecognitionResult ret = new RecognitionResult();
            foreach (Gesture template in temp_templates)
            {
                ret.tests += 1.0f;

                float fscore = template.fscore * template.lower_bound;
                //if (template.lower_bound > template.rejection_threshold)
                if (fscore > template.rejection_threshold)
                {
                    ret.pruned_by_rejection += 1.0f;
                    continue;
                }

                if (fscore > ret.score)
                {
                    ret.pruned_by_best_score += 1.0f;
                    continue;
                }

                float score = DTW_Distance(
                    candidate.vecs,
                    template.vecs,
                    Gesture.r);

                if (score >= template.rejection_threshold)
                    continue;

                score *= template.fscore;
                if (score < ret.score)
                {
                    ret.template = template;
                    ret.score = score;
                }
            }

            return ret;
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

        public static float Clamp(float value, float min, float max)
        {
            return (value < min) ? min : (value > max) ? max : value;
        }
    }
    
    class Gesture
    {        
        
        public static int resample_cnt = 16;
        public static int r = 2;//(resample_cnt / 10);
        
        /*
        public static int resample_cnt = 16;
        public static int r = (resample_cnt / 10);
        */

        public string gname { get; set; }
        public List<Vector<float>> raw_pts;
        public List<Vector<float>> pts;
        public List<Vector<float>> vecs;
        public List<Vector<float>> lower, upper;

        public List<float> rejection;
        public float rejection_threshold;

        internal List<Vector<float>> features;
        internal float lower_bound;
        internal float fscore;

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

            features = Extract_Features(temp);

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

            features = Extract_Features(temp);


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

            features = Extract_Features(temp);

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
            return lower_bound.CompareTo(other.lower_bound);
        }

        static List<Vector<float>> Extract_Features(List<Vector<float>> points)
        {
            if (points.Count == 0)
                return null;

            var m = points[0].Count;
            var abs_dist = Vector<float>.Build.Dense(m, 0);

            //var emin = Vector<float>.Build.Dense(33, float.PositiveInfinity);
            //var emax = Vector<float>.Build.Dense(33, float.NegativeInfinity);

            for (int ii = 1;
                 ii < points.Count;
                 ii++)
            {
                var vec = points[ii] - points[ii - 1];
                Vector<float> point = points[ii];

                for (int jj = 0;
                     jj < m;
                     jj++)
                {
                    //int kk = jj % 33;
                    abs_dist[jj] += Math.Abs(vec[jj]);             
                    //emin[kk] = Math.Min(emin[kk], Math.Abs(vec[jj]));
                    //emax[kk] = Math.Max(emax[kk], Math.Abs(vec[jj]));
                }
            }

            //var edeltas = emax - emin;

            return new List<Vector<float>> {
                abs_dist /  (float) abs_dist.L2Norm(),
                //edeltas / (float) edeltas.L2Norm(),
            };
        }

        public bool Is2D()
        {
            if (gname == "x" || gname == "c" || gname == "circle" || gname == "triangle"
                || gname == "rectangle" || gname == "check" || gname == "caret" || gname == "zigzag"
                || gname == "arrow" || gname == "star")
                return true;
            return false;
        }
    }
}
