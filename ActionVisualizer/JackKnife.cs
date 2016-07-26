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

namespace ActionVisualizer
{
    class JackKnife
    {
        private List<Gesture> templates;

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
        
        public Tuple<Gesture,float> Classify(Gesture candidate)
        {
            Gesture best = null;
            float best_score = float.NegativeInfinity;
            foreach (Gesture template in templates)
            {
                float score = DTW_Distance(candidate, template, (int)(template.vecs.Count * .25f));
                if (score > best_score)
                {
                    best = template;
                    best_score = score;
                }
            }
            return new Tuple<Gesture, float>(best, best_score);
        }

        public float DTW_Distance(Gesture candidate, Gesture template, int w)
        {
            int m = candidate.vecs.Count+1, n = template.vecs.Count+1;
            float[,] dtw = new float[m,n];

            for (int ii = 0; ii < m; ii++)
                for (int jj = 0; jj < n; jj++)
                    dtw[ii, jj] = float.NegativeInfinity;
            dtw[0, 0] = 0;
            for (int ii = 1; ii < m; ii++)
            {
                for (int jj = Math.Max(1, ii - w); jj < Math.Min(m, ii + w); jj++)
                {
                    float cost = candidate.vecs[ii - 1] * template.vecs[jj - 1];
                    float max = Math.Max(dtw[ii - 1, jj], Math.Max(dtw[ii, jj - 1], dtw[ii - 1, jj - 1]));
                    dtw[ii, jj] = cost + max;
                }                   
            }

            return dtw[candidate.vecs.Count, template.vecs.Count]/m;
        }
    }
    
    class Gesture
    {
        static int resample_cnt = 16;
        public string gname { get; set; }
        public List<Vector2> raw_pts;
        public List<Vector2> pts;
        public List<Vector2> vecs;

        public Gesture(List<Vector2> temp, string label)
        {
            gname = label;
            raw_pts = temp;
            pts = Resample(temp,resample_cnt);
            vecs = new List<Vector2>();
            for (int ii = 1; ii < pts.Count; ii++)
            {
                vecs.Add((pts[ii] - pts[ii - 1]).Normalize());
            }
        }

        public Gesture(List<Point> points, string label)
        {
            List<Vector2> temp = new List<Vector2>();
            foreach (Point p in points)
                temp.Add(new Vector2((float)p.X, (float)p.Y));
            gname = label;
            raw_pts = temp;
            pts = Resample(temp, resample_cnt);
            vecs = new List<Vector2>();
            for (int ii = 1; ii < pts.Count; ii++)
            {
                vecs.Add((pts[ii] - pts[ii - 1]).Normalize());
            }
        }

        private float PathLength(List<Vector2> points)
        {
            float length = 0;
            for (int ii = 1; ii < points.Count; ii++)
            {
                length += (points[ii] - points[ii - 1]).Magnitude();
            }
            return length;
        }

        private List<Vector2> Resample(List<Vector2> temp, int n)
        {
            List<Vector2> points = new List<Vector2>();
            foreach (Vector2 v in temp)
                points.Add(new Vector2(v.X, v.Y));
            List<Vector2> ret = new List<Vector2>();
            ret.Add(new Vector2(points[0].X, points[0].Y));
            float I = PathLength(points) / (n - 1.0f);
            float D = 0.0f;
            int ii = 1;
            while (ii < points.Count)
            {
                float d = (points[ii] - points[ii - 1]).Magnitude();
                if (D + d >= I)
                {
                    Vector2 vec = points[ii] - points[ii - 1];
                    float t = (I - D) / d;
                    Vector2 q = points[ii - 1] + vec*t;
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

    }
}
