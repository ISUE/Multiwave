using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Exocortex.DSP;

namespace ActionVisualizer
{
    class KeyFrequency
    {
        public int frequency;
        public int channel;
        public int radius;
        public int prior;
       /*
        public double speakerTheta = 65;
        public double speakerTheta2 = -54;
        */

        public double speakerTheta;
        public double speakerAltitude;
        public float[] data;
        public double[] stackedData;

        public bool isBoth;

        public int state { get; set; }
        public int inverse_state { get; set; }
        public double x;
        public double y;
        public double z;

        public KeyFrequency() { }

        public KeyFrequency(int f, int c, int rad, double[] array) {
            this.frequency = f;
            this.channel = c;
            this.radius = rad;
            this.data = new float[radius * 2 + 1];
            Array.Copy(array, this.data, radius * 2 + 1 );

            isBoth = false;
            inverse_state = 0;
            prior = 0;
        }

        public KeyFrequency(int f, int c, int rad, ComplexF[] indata, int center)
        {
            this.frequency = f;
            this.channel = c;
            this.radius = rad;
            this.data = new float[radius * 2 + 1];

            //rescale numbers >=0

            float min = float.MaxValue;
            float max = float.MinValue;
            for (int i = 0; i < this.data.Length; i++)
            {
                this.data[i] = mag2db(indata[center-rad+i]);
                if (this.data[i] < min)
                    min = this.data[i];
            }
            for (int i = 0; i < this.data.Length; i++)
            {
                this.data[i] -= min;
                if (this.data[i] > max)
                    max = this.data[i];
            }

            max = Math.Max(max, float.Epsilon);

            for (int i = 0; i < this.data.Length; i++)
            {
                this.data[i] /= max;
                if (float.IsNaN(data[i]))
                    Console.WriteLine("NaN");
            }

            isBoth = false;
            inverse_state = 0;
            prior = 0;
            state = this.classify();
        }
        
        public KeyFrequency(int f, int c, int rad, double[] indata, int center, int priori)
        {
            this.frequency = f;
            this.channel = c;
            this.radius = rad;
            this.data = new float[radius * 2 + 1];

            //rescale numbers >=0

            float min = float.MaxValue;
            float max = float.MinValue;
            for (int i = 0; i < this.data.Length; i++)
            {
                this.data[i] = (float)indata[center - rad + i];
                if (this.data[i] < min)
                    min = this.data[i];
            }            

            for (int i = 0; i < this.data.Length; i++)
            {
                this.data[i] -= min;
                if (this.data[i] > max)
                    max = this.data[i];
            }

            max = Math.Max(max, float.Epsilon);
            for (int i = 0; i < this.data.Length; i++)
            {
                this.data[i] /= max;
                if (float.IsNaN(data[i]))
                    Console.WriteLine("NaN");
            }

            isBoth = false;
            inverse_state = 0;
            /*
            if (c == 1)
                speakerTheta = 65;
            else if (c == 2)
                speakerTheta = -54;
            else
                speakerTheta = 0;
              */
            switch (c)
            {
                case 1:
                    speakerTheta = 55;
                    speakerAltitude = -10;
                    break;
                case 2:
                    speakerTheta = -55;
                    speakerAltitude = -10;
                    break;
                case 3:
                    speakerTheta = 0;
                    speakerAltitude = -10;
                    break;
                case 4:
                    speakerTheta = 0;
                    speakerAltitude = 0;
                    break;
                case 5:
                    speakerTheta = 110;
                    speakerAltitude = 30;
                    break;
                case 6:
                    speakerTheta = -110;
                    speakerAltitude = 30;
                    break;

            }
            this.prior = priori;
            state = this.classify();           
        }

        public float mag2db(ComplexF y)
        {
            return 20.0f * (float)Math.Log10(Math.Sqrt(y.Re * y.Re + y.Im * y.Im)/.02);
        }

        public override string ToString()
        {
            string sto = "";
            for (int i = 0; i<this.data.Length; i++)
                sto = sto + " " + (i+1) + ":" + this.data[i];

            return sto;
        }

        //Heuristics classification
        /*public int classify()
        {
            double thresh = .5;
            bool towards = false;
            bool away = false;
            int velocityT = 0;
            int velocityA = 0;
            for (int i = 4; i < radius; i++)
            {
                if (data[radius - i] > thresh)
                {
                    velocityA = i - 3;
                    away = true;
                }
                if (data[radius + i] > thresh)
                {
                    velocityT = i - 3;
                    towards = true;
                }
            }
            if ((towards == true) && (away == true))
            {
                isBoth = true;
                if (velocityA > velocityT)
                {
                    inverse_state = velocityT;
                    updateXY(-1*velocityA);
                    return -1 * velocityA;
                }
                else
                {
                    inverse_state = velocityA;
                    updateXY(velocityT);
                    return velocityT;
                }
            }
            if (towards)
            {
                updateXY(velocityT);
                return velocityT;
            }
            if (away)
            {
                updateXY(-1*velocityA);
                return -1 * velocityA;
            }
            return 0;
        }*/
        public int classify()
        {
            double thresh = .1;
            double second_peak_thresh = .3;

            double peakA = 1.0;
            double peakT = 1.0;

            bool aDropped = true;
            bool tDropped = true;

            bool aSecondPeak = false;
            bool tSecondPeak = false;

            bool towards = false;
            bool away = false;
            int velocityT = 0;
            int velocityA = 0;
            for (int i = 4; i < radius; i++)
            {    
                if (data[radius - i] > thresh)
                {
                    velocityA = i - 3;
                    away = true;
                }
                if (data[radius + i] > thresh)
                {
                    velocityT = i - 3;
                    towards = true;
                }   
            }
            if ((towards == true) && (away == true))
            {
                isBoth = true;
                if (velocityA > velocityT)
                {
                    if (prior > 0)
                    {
                        inverse_state = -1 * velocityA;
                        updateXY(velocityT);
                        return velocityT;
                    }
                    inverse_state = velocityT;
                    updateXY(-1 * velocityA);
                    return -1 * velocityA;
                }
                else
                {
                    if (prior < 0)
                    {
                        inverse_state = velocityT;
                        updateXY(-1 * velocityA);
                        return -1 * velocityA;
                    }
                    inverse_state = -1*velocityA;
                    updateXY(velocityT);
                    return velocityT;
                }
            }
            if (towards)
            {
                updateXY(velocityT);
                return velocityT;
            }
            if (away)
            {
                updateXY(-1 * velocityA);
                return -1 * velocityA;
            }
            return 0;
        }
        private void updateXY(int velocity)
        {
            this.x = (double)velocity * Math.Sin(speakerTheta * Math.PI / 180.0) * Math.Cos(speakerAltitude * Math.PI / 180.0);
            this.y = (double)velocity * Math.Cos(speakerTheta * Math.PI / 180.0) * Math.Cos(speakerAltitude * Math.PI / 180.0);
            this.z = (double)velocity * Math.Sin(speakerAltitude * Math.PI / 180.0);
        }
    }
}
