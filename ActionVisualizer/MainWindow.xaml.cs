﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Documents;
using System.Windows.Ink;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Windows.Shapes;
using Exocortex.DSP;
using NAudio.Dsp;
using NAudio.Wave;
using VerySimpleKalman;



namespace ActionVisualizer
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private AsioOut asioOut;
        private WaveIn waveIn;

        public int waveOutChannels;

        public int minFrequency = 18700;
        public int frequencyStep = 500;
        public int buffersize = 2048;
        public int[] bin;
        public float[] sampledata;
        public float[] inbetween;
        bool init_inbetween = true;
        ComplexF[] indata;

        double[] filteredindata;
        double[] priori;

        int[] channelLabel;
        int[] velocity;
        int[] displacement;

        int[] prev_displacement;
        int[] instant_displacement;
        int[] towards_displacement;

        double ratio;
        VDKalman filter;

        Rectangle[] bars;
        Line[] lines;

        int selectedChannels = 1;
        List<int> frequencies;
        List<int> centerbins;

        List<KeyFrequency> KF;

        List<List<int>> history;
        List<List<int>> inverse_history;
        PointCollection pointHist;
        StylusPointCollection S;

        private readonly List<ModelVisual3D> _models = new List<ModelVisual3D>();
        ModelingHelper modelingHelper;
        Point3DCollection point3DHist;

        ComplexGesture gesture_history;
        ComplexGesture3D gesture_history3D;
        bool readyforgesture = false;
        bool gesture_started = false;
        int motion_free = 0;
        int idle_count = 0;
        int ignoreFrames = 0;

        GestureTests.Logger Log;        

        public MainWindow()
        {
            InitializeComponent();
            this.WindowStartupLocation = System.Windows.WindowStartupLocation.CenterScreen;

            this.KeyDown += new KeyEventHandler(MainWindow_KeyDown);

            int waveInDevices = WaveIn.DeviceCount;
            for (int waveInDevice = 0; waveInDevice < waveInDevices; waveInDevice++)
            {
                WaveInCapabilities deviceInfo = WaveIn.GetCapabilities(waveInDevice);
                Console.WriteLine("Device {0}: {1}, {2} channels",
                    waveInDevice, deviceInfo.ProductName, deviceInfo.Channels);
            }

            waveIn = new WaveIn();
            waveIn.BufferMilliseconds = 47 * buffersize / 2048;
            waveIn.DeviceNumber = 0;
            waveIn.WaveFormat = new WaveFormat(44100, 32, 1);
            waveIn.DataAvailable += waveIn_DataAvailable;

            try
            {
                waveIn.StartRecording();
            }
            catch (NAudio.MmException e)
            {
                Console.WriteLine(e.ToString() + "\nPlug in a microphone!");
            }

            history = new List<List<int>>();
            inverse_history = new List<List<int>>();
            pointHist = new PointCollection();
            point3DHist = new Point3DCollection();
            modelingHelper = new ModelingHelper();

            bin = new int[buffersize * 2];
            sampledata = new float[buffersize * 2];
            priori = new double[buffersize * 2];

            channelLabel = new int[1];
            channelLabel[0] = 1;
            velocity = new int[1];
            velocity[0] = 0;

            prev_displacement = new int[1];
            prev_displacement[0] = 0;

            instant_displacement = new int[1];
            instant_displacement[0] = 0;

            towards_displacement = new int[1];
            towards_displacement[0] = 1;


            displacement = new int[1];
            displacement[0] = 0;
            bars = new Rectangle[1];
            bars[0] = createBar(channelLabel[0], selectedChannels, velocity[0], 33);
            lines = new Line[1];
            lines[0] = createLine(channelLabel[0], selectedChannels);
            for (int i = 0; i < buffersize * 2; i++)
            {
                bin[i] = i;
                sampledata[i] = 0;
                priori[i] = 0;

            }

            _barcanvas.Background = new SolidColorBrush(Colors.Black);

            filter = new VDKalman(2);
            filter.initialize(1, .1, 1, 0);

            history.Add(new List<int> { 0 });
            inverse_history.Add(new List<int> { 0 });

            Log = new GestureTests.Logger("ActionVisualizer");
            gesture_history = new ComplexGesture();
            gesture_history3D = new ComplexGesture3D();
            WekaHelper.initialize();
        }

        void waveIn_DataAvailable(object sender, WaveInEventArgs e)
        {
            //Console.WriteLine("WaveIn_DataAvailable");
            //Console.WriteLine(e.BytesRecorded); //8288 bytes -> 2072 floats (24 too many)
            if (init_inbetween)
            {
                inbetween = new float[e.BytesRecorded / 4 - buffersize];
                for (int i = 0; i < e.BytesRecorded / 4 - buffersize; i++)
                    inbetween[i] = 0;
                init_inbetween = false;
            }


            for (int index = 0; index < buffersize; index++)
            {
                int sample = (int)((e.Buffer[index * 4 + 3] << 24) |
                                        (e.Buffer[index * 4 + 2] << 16) |
                                        (e.Buffer[index * 4 + 1] << 8) |
                                        (e.Buffer[index * 4 + 0]));

                float sample32 = sample / 2147483648f;

                if (index >= (buffersize - inbetween.Length))
                    sampledata[index] = inbetween[index - buffersize + inbetween.Length];
                else
                    sampledata[index] = sampledata[index + buffersize + inbetween.Length];
                sampledata[index + buffersize] = sample32;

            }

            if (e.BytesRecorded / 4 - buffersize < 0)
                return;
            inbetween = new float[e.BytesRecorded / 4 - buffersize];

            for (int i = buffersize; i < e.BytesRecorded / 4; i++)
            {
                int sample = (int)((e.Buffer[i * 4 + 3] << 24) |
                                        (e.Buffer[i * 4 + 2] << 16) |
                                        (e.Buffer[i * 4 + 1] << 8) |
                                        (e.Buffer[i * 4 + 0]));

                float sample32 = sample / 2147483648f;
                inbetween[i - buffersize] = sample32;
            }
            /*
            for (int i = 0; i < 4096; ++i)
                sampledata[i] = 0.01*(float)Math.Sin(2 * Math.PI * i * 19000 / 44100f);
            */
            bufferFFT();
            filter.time_Update();
            if ((asioOut != null))
            {
                foreach (Rectangle r in bars)
                    _barcanvas.Children.Remove(r);

                KF = new List<KeyFrequency>();
                //gestureDetected.Text = "";
                for (int i = 0; i < frequencies.Count; i++)
                {
                    if (history[i].Count > 0)
                        KF.Add(new KeyFrequency(frequencies.ElementAt(i), i + 1, frequencyStep/30, filteredindata, centerbins.ElementAt(i), history[i].Last()));
                    else
                        KF.Add(new KeyFrequency(frequencies.ElementAt(i), i + 1, frequencyStep/30, filteredindata, centerbins.ElementAt(i), 0));

                    velocity[i] = KF.ElementAt(i).state;
                    filter.measurement_Update(velocity[i]);
                    displacement[i] += (int)filter.x_priori[0];

                    //gestureDetected.Text += " " + displacement[i];
                    bars[i] = createBar(channelLabel[i], selectedChannels, velocity[i], frequencyStep/30+3);
                    _barcanvas.Children.Add(bars[i]);

                    if ((displacement[i] - prev_displacement[i]) < 0 && velocity[i] > 0)
                    {
                        //Console.WriteLine("Away Stop: " + (instant_displacement[i]-prev_displacement[i]));
                        ratio = Math.Abs((double)((instant_displacement[i] - prev_displacement[i]) / (double)towards_displacement[i]));
                        //double ratio = (double)((instant_displacement[i] - prev_displacement[i]) + (double)towards_displacement[i]);

                        //Console.WriteLine(filter.x_priori[0]);
                        prev_displacement[i] = displacement[i];
                    }
                    if ((displacement[i] - prev_displacement[i]) > 0 && velocity[i] < 0)
                    {
                        //Console.WriteLine("Towards Stop: " + (instant_displacement[i] - prev_displacement[i]));
                        towards_displacement[i] = instant_displacement[i] - prev_displacement[i];
                        prev_displacement[i] = displacement[i];
                    }
                    instant_displacement[i] = displacement[i];

                    history[i].Add(velocity[i]);
                    inverse_history[i].Add(KF[i].inverse_state);
                    //if (history[i].Count > 20)
                    //    history[i].RemoveAt(0);
                }
                detectGestures();
            }
        }

        private void detectGestures()
        {
            ignoreFrames++;
            if (selectedChannels == 1)
            {
                foreach (List<int> subList in history)
                {
                    int signChanges = 0, bandwidth = 0, step = 0, lastSig = 0;
                    for (int i = 1; i < subList.Count; i++)
                    {
                        step++;
                        if (subList[i - 1] != 0)
                            lastSig = subList[i - 1];
                        if (subList[i] * lastSig < 0)
                        {
                            signChanges++;
                            bandwidth += step;
                            step = 0;
                        }
                    }
                    
                    if (KF[0].isBoth && KF[0].inverse_state > 5)
                        gestureDetected.Text = "Two Handed ";
                    else if (signChanges == 0 && (lastSig != 0))
                        gestureDetected.Text = "Scrolling ";
                    else if (signChanges == 2 || signChanges == 1)
                        gestureDetected.Text = "SingleTap ";
                    else if (signChanges >= 3)
                        gestureDetected.Text = "DoubleTap ";
                }
            }
            else if (selectedChannels == 2)
            {
                double tot_X = 0, tot_Y = 0;
                foreach (KeyFrequency now in KF)
                {
                    tot_X += now.x;
                    tot_Y += now.y;
                }
                //gestureDetected.Text += (Math.Round(tot_X,2) + " " + Math.Round(tot_Y,2));
                pointHist.Add(new Point(tot_X, tot_Y));
                //if (pointHist.Count > 20)
                //    pointHist.RemoveAt(0);
                if (!gesture_started && tot_X == 0 && tot_Y == 0)
                {
                    idle_count++;

                    pointHist.Clear();
                    foreach (List<int> sublist in history)
                        sublist.Clear();
                    foreach (List<int> sublist in inverse_history)
                        sublist.Clear();


                    //TODO, Evaluate the proper length of the deadzone. (seems to be 6 frames)
                    if (idle_count > frameWindow.SelectedIndex)
                    {
                        idle_count = 0;
                        gesture_history.clearGestureHistory();
                        //gestureDetected.Text = "none";
                    }
                }
                if (gesture_started && tot_X == 0 && tot_Y == 0)
                    motion_free++;
                if (tot_X != 0 || tot_Y != 0)
                {
                    gesture_started = true;
                    motion_free = 0;
                    idle_count = 0;
                }

                generateStroke(pointHist);
            }
            else if (selectedChannels >= 3)
            {
                double tot_X = 0, tot_Y = 0, tot_Z = 0;
                foreach (KeyFrequency now in KF)
                {

                    tot_X += now.x;
                    tot_Y += now.y;
                    tot_Z += now.z;
                    //For Graphing
                    //tot_X -= now.x;
                    //tot_Y += now.z;
                    //tot_Z -= now.y;
                }
                //gestureDetected.Text += (Math.Round(tot_X, 2) + " " + Math.Round(tot_Y, 2) + " " + Math.Round(tot_Z, 2));
                pointHist.Add(new Point(tot_X, tot_Y));
                point3DHist.Add(new Point3D(tot_X, tot_Y, tot_Z));
                //if (pointHist.Count > 20)
                //    pointHist.RemoveAt(0);
                //if (point3DHist.Count > 20)
                //    point3DHist.RemoveAt(0);
                if (!gesture_started && tot_X == 0 && tot_Y == 0 && tot_Z == 0)
                {
                    idle_count++;
                    pointHist.Clear();
                    point3DHist.Clear();
                    foreach (List<int> sublist in history)
                        sublist.Clear();
                    foreach (List<int> sublist in inverse_history)
                        sublist.Clear();

                    //TODO, Evaluate the proper length of the deadzone. 
                    if (idle_count > frameWindow.SelectedIndex)
                    {
                        idle_count = 0;
                        gesture_history3D.clearGestureHistory();
                        //gestureDetected.Text = "none";
                    }
                }
                if (gesture_started && tot_X == 0 && tot_Y == 0 && tot_Z == 0)
                    motion_free++;
                if (tot_X != 0 || tot_Y != 0 || tot_Z != 0)
                {
                    gesture_started = true;
                    motion_free = 0;
                    idle_count = 0;
                }

                generateStroke(pointHist);
                //addPoints(point3DHist);
            }
            
            
            gestureCompleted();
        }

        public void generateStroke(PointCollection pointHist)
        {
            S = new StylusPointCollection();
            S.Add(new StylusPoint(_ink.ActualWidth / 2, _ink.ActualHeight / 2));
            for (int i = 0; i < pointHist.Count; i++)
            {
                S.Add(new StylusPoint(S[i].X - pointHist[i].X, S[i].Y - pointHist[i].Y));
            }
            Stroke So = new Stroke(S);
            _ink.Strokes.Clear();
            _ink.Strokes.Add(So);

        }

        public void bufferFFT()
        {
            indata = new ComplexF[buffersize * 2];
            for (int i = 0; i < buffersize * 2; i++)
            {
                indata[i].Re = sampledata[i] * (float)FastFourierTransform.HammingWindow(i, buffersize * 2);
                indata[i].Im = 0;
            }
            Exocortex.DSP.Fourier.FFT(indata, buffersize * 2, Exocortex.DSP.FourierDirection.Forward);
            filteredindata = filterMean(indata, 1.3);
        }

        private void StartStopSineWave()
        {
            if (asioOut == null)
            {
                button1.Content = "Stop Sound";
                Console.WriteLine("User Selected Channels: " + selectedChannels);
                WaveOutCapabilities outdeviceInfo = WaveOut.GetCapabilities(0);
                waveOutChannels = outdeviceInfo.Channels;
                asioOut = new AsioOut(0);
                Console.WriteLine("Asio Driver Output Count: " + asioOut.DriverName);
                
                int waveOutDevices = WaveOut.DeviceCount;
                for (int i = 0; i < waveOutDevices; i++)
                {
                    outdeviceInfo = WaveOut.GetCapabilities(i);
                    Console.WriteLine("Device {0}: {1}, {2} channels",
                            i, outdeviceInfo.ProductName, outdeviceInfo.Channels);                    
                }

                List<IWaveProvider> inputs = new List<IWaveProvider>();
                frequencies = new List<int>();
                centerbins = new List<int>();

                for (int c = 0; c < selectedChannels; c++)
                {
                    if (c >=4)
                    {
                        //Make LR/RR speakers fit in frequency space.
                        inputs.Add(new SineWaveProvider32(minFrequency + ((c - 1) * frequencyStep), 1.0f, 44100, 1));
                        //inputs.Add(new SineWaveProvider32(18000 + c * 700, 1f, 44100, 1));
                        frequencies.Add(minFrequency + ((c - 1) * frequencyStep));
                        centerbins.Add((int)Math.Round((minFrequency + ((c - 1) * frequencyStep)) / 10.768));    
                    }
                    else if (c >= 3)
                    {
                        //Make LR/RR speakers fit in frequency space.
                        inputs.Add(new SineWaveProvider32(minFrequency + ((c - 1) * frequencyStep), 0.0f, 44100, 1));
                        //inputs.Add(new SineWaveProvider32(18000 + c * 700, 1f, 44100, 1));
                        frequencies.Add(minFrequency + ((c - 1) * frequencyStep));
                        centerbins.Add((int)Math.Round((minFrequency + ((c - 1) * frequencyStep)) / 10.768));                        
                    }
                    else
                    {
                        //Original Sine Wave generation
                        inputs.Add(new SineWaveProvider32(minFrequency + c * frequencyStep, 0.5f, 44100, 1));
                        frequencies.Add(minFrequency + c * frequencyStep);
                        centerbins.Add((int)Math.Round((minFrequency + c * frequencyStep) / 10.768));
                    }
                }

                var splitter = new MultiplexingWaveProvider(inputs, selectedChannels);
                try
                {
                    asioOut.Init(splitter);
                    asioOut.Play();
                }
                catch (System.ArgumentException)
                {
                    Console.WriteLine("Invalid audio channel count. Please select a lower number of audio channels");
                }
                //waveOut = new WaveOut();
                //waveOut.Init(sineWaveProvider);                    
                //waveOut.Init(splitter);

                Console.WriteLine("Number of Channels: " + asioOut.NumberOfOutputChannels);
                Console.WriteLine("Playback Latency: " + asioOut.PlaybackLatency);
            }
            else
            {
                asioOut.Stop();
                asioOut.Dispose();
                asioOut = null;
                button1.Content = "Start Sound";

                frequencies.Clear();
                centerbins.Clear();
                //Foutstream.Clear();
            }
        }

        public float mag2db(ComplexF y)
        {
            return 20.0f * (float)Math.Log10(Math.Sqrt(y.Re * y.Re + y.Im * y.Im) / .02);
        }

        public double calculateTotalAngleXY(StylusPointCollection points)
        {
            double TotalAngleXY = 0.0;

            int N = points.Count;

            double dx, dy;

            for (int i = 1; i < N; ++i)
            {
                dx = points[i].X - points[i - 1].X;
                dy = points[i].Y - points[i - 1].Y;

                float angleXY = (float)Math.Atan2(dy, dx);

                TotalAngleXY += angleXY*angleXY;
            }

            return TotalAngleXY;
        }

        public double[] filterMean(ComplexF[] data, double factor)
        {
            double[] outdata = new double[data.Length];
            if (false)//centerbins != null)
            {
                for (int c = 0; c < centerbins.Count; c++)
                {
                    double min = Double.PositiveInfinity;
                    double mean = 0;
                    for (int i = centerbins[c] - 33; i <= (centerbins[c] + 33); i++)
                    {
                        outdata[i] = mag2db(data[i]);
                        min = Math.Min(outdata[i], min);
                    }

                    for (int i = centerbins[c] - 33; i <= (centerbins[c] + 33); i++)
                    {
                        //outdata[i] -= min;
                        mean += (outdata[i]);
                    }
                }
            }
            else
            {
                double min = Double.PositiveInfinity;
                double mean = 0;
                for (int i = 0; i < data.Length; i++)
                {
                    outdata[i] = mag2db(data[i]);
                    min = Math.Min(outdata[i], min);
                }

                for (int i = 0; i < data.Length; i++)
                {
                    outdata[i] -= min;
                    mean += (outdata[i]);
                }
                mean /= data.Length;
                for (int i = 0; i < data.Length; i++)
                    if (outdata[i] < (mean * factor))
                        outdata[i] = 0;

                for (int i = 0; i < data.Length; i++)
                {
                    if ((i > 0) && (i < (data.Length - 1)))
                    {
                        if ((outdata[i] > 0) && (priori[i] == 0) && (outdata[i - 1] == 0) && (outdata[i + 1] == 0))
                        {
                            outdata[i] = 0;
                        }
                    }
                    priori[i] = outdata[i];
                }
            }
            return outdata;
        }

        private void channelSelector_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (asioOut != null)
                StartStopSineWave();

            selectedChannels = (sender as ComboBox).SelectedIndex + 1;

            if (selectedChannels == 2)
            {
                _ink.Visibility = Visibility.Visible;
                _viewport.Visibility = Visibility.Hidden;
            }
            if (selectedChannels == 6)
            {
                _ink.Visibility = Visibility.Visible;
                //_viewport.Visibility = Visibility.Visible;
            }

            history.Clear();
            inverse_history.Clear();

            channelLabel = new int[selectedChannels];
            velocity = new int[selectedChannels];
            displacement = new int[selectedChannels];

            foreach (Rectangle l in bars)
                _barcanvas.Children.Remove(l);

            bars = new Rectangle[selectedChannels];

            instant_displacement = new int[selectedChannels];
            prev_displacement = new int[selectedChannels];
            towards_displacement = new int[selectedChannels];

            foreach (Line l in lines)
                _barcanvas.Children.Remove(l);

            lines = new Line[selectedChannels + 1];

            for (int i = 0; i < selectedChannels; i++)
            {
                history.Add(new List<int> { 0 });
                inverse_history.Add(new List<int> { 0 });
                channelLabel[i] = i + 1;
                velocity[i] = 0;

                instant_displacement[i] = 0;
                prev_displacement[i] = 0;
                towards_displacement[i] = 1;

                displacement[i] = 0;
                bars[i] = createBar(i + 1, selectedChannels, 0, 33);
            }
            for (int i = 0; i < lines.Length; i++)
            {
                lines[i] = createLine(i + 1, selectedChannels);
                _barcanvas.Children.Add(lines[i]);
            }

        }

        private Rectangle createBar(int channel, int cmax, int velocity, int vmax)
        {
            Rectangle outRect = new Rectangle();
            outRect.Width = _barcanvas.ActualWidth / (cmax + 2);
            outRect.Height = (_barcanvas.ActualHeight * Math.Abs(velocity)) / (2 * vmax);
            outRect.Stroke = new SolidColorBrush(Colors.Black);
            Canvas.SetLeft(outRect, (outRect.Width * (channel)));

            if (velocity < 0)
            {
                Canvas.SetTop(outRect, _barcanvas.ActualHeight / 2);
                outRect.Fill = new SolidColorBrush(Colors.Red);
            }
            else
            {
                Canvas.SetTop(outRect, (_barcanvas.ActualHeight / 2) - ((_barcanvas.ActualHeight * Math.Abs(velocity)) / (2 * vmax)));
                outRect.Fill = new SolidColorBrush(Colors.Blue);
            }

            //Console.WriteLine(outRect.Width + " " + outRect.Height);
            return outRect;
        }

        private Line createLine(int channel, int cmax)
        {
            double offset = _barcanvas.ActualWidth / (cmax + 2);
            Line outLine = new Line()
            {
                X1 = _barcanvas.ActualWidth * channel / (cmax + 2),
                Y1 = 0,
                X2 = _barcanvas.ActualWidth * channel / (cmax + 2),
                Y2 = _barcanvas.ActualHeight,
                StrokeThickness = 1,
                Stroke = new SolidColorBrush(Colors.White)
                //StrokeDashArray = 
            };
            //            Canvas.SetLeft(outLine, offset * channel);
            return outLine;
        }

        private void button1_Click(object sender, RoutedEventArgs e)
        {
            StartStopSineWave();
        }

        public void addPoints(Point3DCollection p)
        {
            _models.Clear();
            _viewport.Children.Clear();
            foreach (Point3D Ps in p)
                _models.Add(modelingHelper.CreateSphere(Ps, 1, 5, 5, Colors.Green));
            _models.ForEach(x => _viewport.Children.Add(x));
        }

        public void gestureCompleted()
        {
            int motion_threshold = 3; //originally 5         
            int ignore_threshold = 20;
             
            if( ignoreFrames <= ignore_threshold)                          
            {
                motion_free = 0;
                readyforgesture = false;
                colorBox.Background = new SolidColorBrush(Colors.Red);
                gesture_started = false;
                //Clear the buffers
                foreach (List<int> sublist in history)
                    sublist.Clear();
                foreach (List<int> sublist in inverse_history)
                    sublist.Clear();
                pointHist.Clear();
                point3DHist.Clear();
            }

            if (gesture_started && ignoreFrames > ignore_threshold && motion_free > motion_threshold && selectedChannels >= 2)
            {

                pointHist = new PointCollection(pointHist.Reverse().Skip(motion_threshold).Reverse());
                point3DHist = new Point3DCollection(point3DHist.Reverse().Skip(motion_threshold).Reverse());
                S = new StylusPointCollection(S.Reverse().Skip(motion_threshold).Reverse());               
                for (int i = 0; i < history.Count; i++)
                {
                    history[i] = new List<int>(history[i].Reverse<int>().Skip(motion_threshold).Reverse<int>());
                    inverse_history[i] = new List<int>(inverse_history[i].Reverse<int>().Skip(motion_threshold).Reverse<int>());
                }

                if (/*gestureSelector.Text == "Detect"*/ detectMode.IsChecked.Value && pointHist.Count > 2)
                {
                    //Call function to find features and test with weka machine
                    if (selectedChannels == 2)
                    {
                        float[] speakers = { (float)KF[0].speakerTheta, (float)KF[1].speakerTheta };
                        string temp = WekaHelper.Classify(useRubine.IsChecked.Value, pointHist.Count() * waveIn.BufferMilliseconds,
                            handSelector.Text == "right", new List<float>(speakers), pointHist, S, history, inverse_history);
                        switch(temp)
                        {
                            case "swipe_up":
                                gestureDetected.Text = "swipe_forward";
                                break;
                            case "swipe_down":
                                gestureDetected.Text = "swipe_back";
                                break;
                            case "tap_up":
                                gestureDetected.Text = "tap_forward";
                                break;
                            case "tap_down":
                                gestureDetected.Text = "tap_back";
                                break;                            
                            default:
                                gestureDetected.Text = temp;
                                break;
                        }                        

                    }
                    if (selectedChannels == 3)
                    {
                        float[] speakers = { (float)KF[0].speakerTheta, (float)KF[1].speakerTheta, (float)KF[2].speakerTheta };
                        float[] elevations = { (float)KF[0].speakerAltitude, (float)KF[1].speakerAltitude, (float)KF[2].speakerAltitude };
                        gestureDetected.Text = WekaHelper.Classify3D(useRubine.IsChecked.Value, pointHist.Count() * waveIn.BufferMilliseconds,
                            handSelector.Text == "right", new List<float>(speakers), new List<float>(elevations), point3DHist, history, inverse_history);
                    }
                    if (selectedChannels == 6)
                    {
                        float[] speakers = { (float)KF[0].speakerTheta, (float)KF[1].speakerTheta, (float)KF[2].speakerTheta, (float)KF[3].speakerTheta, (float)KF[4].speakerTheta, (float)KF[5].speakerTheta };
                        float[] elevations = { (float)KF[0].speakerAltitude, (float)KF[1].speakerAltitude, (float)KF[2].speakerAltitude, (float)KF[3].speakerAltitude, (float)KF[4].speakerAltitude, (float)KF[5].speakerAltitude };
                        gestureDetected.Text = WekaHelper.Classify3D(useRubine.IsChecked.Value, pointHist.Count() * waveIn.BufferMilliseconds,
                            handSelector.Text == "right", new List<float>(speakers), new List<float>(elevations), point3DHist, history, inverse_history);
                    }

                    //All the parameters to be passed to ComplexGesture are passed here.
                    double deltaX = S.Last().X - S.First().X;
                    double deltaY = S.Last().Y - S.First().Y;

                    double magnitude = Math.Sqrt(deltaX * deltaX + deltaY * deltaY);
                    double angle = Math.Atan2(deltaY,deltaX);

                    double duration = pointHist.Count() * waveIn.BufferMilliseconds;

                    //Add gestures to the gesture history.
                    if (selectedChannels == 2)
                    {
                        gesture_history.addGesture(angle, magnitude, duration, calculateTotalAngleXY(S), gestureDetected.Text);
                        string gest = gesture_history.updateDetectedGesture();
                        complexGestureDetected.Text = gest;
                        if(ComplexGesture.isComplex(gest))
                        {
                            ignoreFrames = 0;
                        }

                    }
                    else if (selectedChannels == 6)
                    {
                        double deltaZ = (from p3d in point3DHist select p3d.Z).Sum();
                        gesture_history3D.addGesture(angle, magnitude, Math.Atan2(deltaZ, magnitude) /*elevation*/, duration, calculateTotalAngleXY(S), gestureDetected.Text);
                        string gest3D = gesture_history3D.updateDetectedGesture();
                        complexGestureDetected.Text = gest3D;
                        if(ComplexGesture3D.isComplex(gest3D))
                        {
                            ignoreFrames = 0;
                        }

                    }

                    if(gestureDetected.Text == gestureSelector.Text)
                        Log.Log(gestureDetected.Text, gestureSelector.Text, history);                    

                    Log.Log(complexGestureDetected.Text, gestureSelector.Text, history);                    
                }
                else if (readyforgesture && pointHist.Count > 2)
                {
                    if (selectedChannels == 2)
                        writeTo2DFile();
                    else
                        writeTo3DFile();
                }
                //Clear the buffers
                foreach (List<int> sublist in history)
                    sublist.Clear();
                foreach (List<int> sublist in inverse_history)
                    sublist.Clear();
                pointHist.Clear();
                point3DHist.Clear();

                //prepare for next gesture (might need a button press)
                readyforgesture = false;
                colorBox.Background = new SolidColorBrush(Colors.Red);
                gesture_started = false;
                motion_free = 0;
            }
        }

        public void writeTo2DFile()
        {
            string DataPath = @"..\..\..\data\u010\";
            string ImagePath = @"..\..\..\image\u001\";
            string StrokePath = @"..\..\..\stroke\u001\";
            string searchPattern = gestureSelector.Text + "???";
            DirectoryInfo di = new DirectoryInfo(DataPath);
            FileInfo[] files = di.GetFiles(searchPattern);
            int file_index = files.Length;

            //CreateSaveBitmap(_ink, ImagePath + gestureSelector.Text + file_index + ".png");
            //SaveStroke(S, StrokePath + gestureSelector.Text + file_index + ".isf");

            //CreateSaveGroupStroke(_ink, ImagePath + gestureSelector.Text + file_index + "_group.png", StrokePath + gestureSelector.Text + file_index + "_group.isf");

            string filename = DataPath + gestureSelector.Text + file_index;
            StreamWriter file = File.CreateText(filename);

            file.WriteLine("GestureName: " + gestureSelector.Text);
            file.WriteLine("Duration(ms): " + pointHist.Count() * waveIn.BufferMilliseconds);
            file.WriteLine("Handedness: " + handSelector.Text);
            file.WriteLine();
            file.WriteLine("SpeakerAngles: " + selectedChannels);
            for (int i = 0; i < KF.Count; i++)
                file.WriteLine(KF[i].speakerTheta);

            file.WriteLine();
            file.WriteLine("InterpretedPoints: " + pointHist.Count());
            foreach (Point p in pointHist)
                file.WriteLine(p.X + "," + p.Y);

            file.WriteLine();
            file.WriteLine("StrokePoints: " + S.Count());
            foreach (StylusPoint p in S)
                file.WriteLine(p.X + "," + p.Y);

            file.WriteLine();
            file.WriteLine("Velocities: " + history[0].Count);
            for (int i = 0; i < history[0].Count; i++)
                file.WriteLine(history[0][i] + "," + history[1][i]);

            file.WriteLine();
            file.WriteLine("InverseVelocities: " + inverse_history[0].Count);
            for (int i = 0; i < inverse_history[0].Count; i++)
                file.WriteLine(inverse_history[0][i] + "," + inverse_history[1][i]);

            file.Flush();
            file.Close();
        }

        public void writeTo3DFile()
        {
            string DataPath = @"..\..\..\data6D\u010\";
            string searchPattern = gestureSelector.Text + "???";
            DirectoryInfo di = new DirectoryInfo(DataPath);
            FileInfo[] files = di.GetFiles(searchPattern);
            int file_index = files.Length;

            string filename = DataPath + gestureSelector.Text + file_index;
            gestureDetected.Text = gestureSelector.Text+ file_index;
            StreamWriter file = File.CreateText(filename);

            file.WriteLine("GestureName: " + gestureSelector.Text);
            file.WriteLine("Duration(ms): " + pointHist.Count() * waveIn.BufferMilliseconds);
            file.WriteLine("Handedness: " + handSelector.Text);
            file.WriteLine();
            file.WriteLine("SpeakerAngles: " + selectedChannels);
            for (int i = 0; i < KF.Count; i++)
                file.WriteLine(KF[i].speakerTheta);
            file.WriteLine();

            file.WriteLine("SpeakerElevations: " + selectedChannels);
            for (int i = 0; i < KF.Count; i++)
                file.WriteLine(KF[i].speakerAltitude);

            file.WriteLine();
            file.WriteLine("InterpretedPoints: " + pointHist.Count());
            foreach (Point3D p in point3DHist)
                file.WriteLine(p.X + "," + p.Y + "," + p.Z);

            file.WriteLine();
            file.WriteLine("StrokePoints: " + S.Count());
            Point3D origin = new Point3D(0, 0, 0);
            file.WriteLine(origin.X + "," + origin.Y + "," + origin.Z);
            for (int i = 0; i < point3DHist.Count(); i++)
            {
                origin.X += point3DHist[i].X;
                origin.Y += point3DHist[i].Y;
                origin.Z += point3DHist[i].Z;
                file.WriteLine(origin.X + "," + origin.Y + "," + origin.Z);
            }

            file.WriteLine();
            file.WriteLine("Velocities: " + history[0].Count);
            for (int i = 0; i < history[0].Count; i++)
                file.WriteLine(history[0][i] + "," + history[1][i] + "," + history[2][i]);

            file.WriteLine();
            file.WriteLine("InverseVelocities: " + inverse_history[0].Count);
            for (int i = 0; i < inverse_history[0].Count; i++)
                file.WriteLine(inverse_history[0][i] + "," + inverse_history[1][i] + "," + inverse_history[2][i]);

            file.Flush();
            file.Close();
        }

        void MainWindow_KeyDown(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.A)
            {
                readyforgesture = true;
                colorBox.Background = new SolidColorBrush(Colors.Green);
            }
            if (e.Key == Key.C)
            {
                _ink.Strokes.Clear();
            }

            if (e.Key == Key.M)
            {
                generateAllGestureStrokes();
            }
        }

        private void SaveStroke(StylusPointCollection S, string filename)
        {
            StrokeCollection col = new StrokeCollection();
            col.Add(new Stroke(S));

            FileStream file = File.Create(filename);
            col.Save(file);
            file.Flush();
            file.Close();
        }

        private void CreateSaveBitmap(InkCanvas canvas, string filename)
        {
            Rect bounds = VisualTreeHelper.GetDescendantBounds(canvas);
            double dpi = 96d;


            RenderTargetBitmap rtb = new RenderTargetBitmap((int)bounds.Width, (int)bounds.Height, dpi, dpi, System.Windows.Media.PixelFormats.Default);


            DrawingVisual dv = new DrawingVisual();
            using (DrawingContext dc = dv.RenderOpen())
            {
                VisualBrush vb = new VisualBrush(canvas);
                dc.DrawRectangle(vb, null, new Rect(new Point(), bounds.Size));
            }

            rtb.Render(dv);

            BitmapEncoder pngEncoder = new PngBitmapEncoder();
            pngEncoder.Frames.Add(BitmapFrame.Create(rtb));

            try
            {
                System.IO.MemoryStream ms = new System.IO.MemoryStream();

                pngEncoder.Save(ms);
                ms.Close();

                System.IO.File.WriteAllBytes(filename, ms.ToArray());
            }
            catch (Exception err)
            {
                MessageBox.Show(err.ToString(), "Error", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        private void CreateSaveGroupStroke(InkCanvas _ink, string img_file, string str_file)
        {
            StrokeCollection InkGroup = new StrokeCollection();
            for (int j = 0; j < pointHist.Count; j++)
            {
                StylusPointCollection T = new StylusPointCollection();
                T.Add(new StylusPoint(_ink.ActualWidth / 2, _ink.ActualHeight / 2));
                for (int i = j; i < pointHist.Count; i++)
                {
                    T.Add(new StylusPoint(T[i - j].X - pointHist[i].X, T[i - j].Y - pointHist[i].Y));
                }
                Stroke So = new Stroke(T);
                InkGroup.Add(So);
            }

            FileStream file = File.Create(str_file);
            InkGroup.Save(file);
            file.Flush();
            file.Close();

            _ink.Strokes.Clear();
            _ink.Strokes.Add(InkGroup);
            CreateSaveBitmap(_ink, img_file);
        }

        private void generateAllGestureStrokes()
        {
            string ImagePath = @"..\..\..\image\u001\";
            string StrokePath = @"..\..\..\stroke\u001\";

            for (int i = 0; i < gestureSelector.Items.Count; i++)
            {
                string searchPattern = ((ComboBoxItem)gestureSelector.Items[i]).Content + "???.isf";
                string GroupSearchPattern = ((ComboBoxItem)gestureSelector.Items[i]).Content + "*";
                DirectoryInfo di = new DirectoryInfo(StrokePath);
                FileInfo[] files = di.GetFiles(searchPattern);
                FileInfo[] group_files = di.GetFiles(GroupSearchPattern);

                if (files.Length == 0)
                    continue;

                _ink.Strokes.Clear();
                StrokeCollection inStrokes = new StrokeCollection();
                foreach (FileInfo file in files)
                {
                    FileStream temp = file.Open(FileMode.Open);
                    inStrokes.Add(new StrokeCollection(temp));
                    temp.Close();
                }
                _ink.Strokes.Add(inStrokes);

                CreateSaveBitmap(_ink, ImagePath + ((ComboBoxItem)gestureSelector.Items[i]).Content + "_all.png");

                _ink.Strokes.Clear();
                inStrokes = new StrokeCollection();
                foreach (FileInfo file in group_files)
                {
                    FileStream temp = file.Open(FileMode.Open);
                    inStrokes.Add(new StrokeCollection(temp));
                    temp.Close();
                }
                _ink.Strokes.Add(inStrokes);
                CreateSaveBitmap(_ink, ImagePath + ((ComboBoxItem)gestureSelector.Items[i]).Content + "_all_group.png");

            }
            _ink.Strokes.Clear();
            generateHTML();
        }

        private void generateHTML()
        {
            string ImagePath = @"..\..\..\image\u001\";

            StreamWriter outfile = new StreamWriter(@"..\..\..\image\u001\visual.html");
            outfile.Write("<!DOCTYPE html><html><body><table style=\"width:300px\"><tr><th></th>" +
                "<th>1</th><th>2</th><th>3</th><th>4</th><th>5</th><th>6</th>" +
                "<th>7</th><th>8</th><th>9</th><th>10</th><th>11</th><th>12</th>" +
                "<th>13</th><th>14</th><th>15</th><th>16</th><th>17</th><th>18</th>" +
                "<th>19</th><th>20</th><th>21</th><th>22</th><th>23</th><th>24</th><th>25</th></tr>");

            for (int i = 0; i < gestureSelector.Items.Count; i++)
            {
                string searchPattern = ((ComboBoxItem)gestureSelector.Items[i]).Content + "?.png";
                string GroupSearchPattern = ((ComboBoxItem)gestureSelector.Items[i]).Content + "?_group.png";

                if (i == 12)
                    GroupSearchPattern = "c?_group.png";

                DirectoryInfo di = new DirectoryInfo(ImagePath);
                FileInfo[] files = di.GetFiles(searchPattern);
                FileInfo[] group_files = di.GetFiles(GroupSearchPattern);
                Console.WriteLine(files.Length);
                Console.WriteLine(group_files.Length);
                if (files.Length == 0)
                    continue;

                outfile.Write("<tr><th>" + ((ComboBoxItem)gestureSelector.Items[i]).Content + "</th>");

                foreach (FileInfo file in files)
                {
                    outfile.Write("<th>" + "<img src=\"" + file.Name + "\">" + "</th>");
                }

                outfile.Write("<th>" + "<img src=\"" + ((ComboBoxItem)gestureSelector.Items[i]).Content + "_all.png\">" + "</th></tr>");

                outfile.Write("<tr><th>" + ((ComboBoxItem)gestureSelector.Items[i]).Content + "_group</th>");

                foreach (FileInfo file in group_files)
                {
                    outfile.Write("<th>" + "<img src=\"" + file.Name + "\">" + "</th>");
                }
                outfile.Write("<th>" + "<img src=\"" + ((ComboBoxItem)gestureSelector.Items[i]).Content + "_all_group.png\">" + "</th></tr>");
            }
            outfile.Write("</table></body></html>");
            outfile.Close();
        }

        private void use3DGestures_Checked(object sender, RoutedEventArgs e)
        {
            gestureSelector.Items.Clear();
            gestureSelector.Items.Add("swipe_left");
            gestureSelector.Items.Add("swipe_right");
            gestureSelector.Items.Add("swipe_up");
            gestureSelector.Items.Add("swipe_down");
            gestureSelector.Items.Add("swipe_front");
            gestureSelector.Items.Add("swipe_back");
            gestureSelector.Items.Add("tap_left");
            gestureSelector.Items.Add("tap_right");
            gestureSelector.Items.Add("tap_up");
            gestureSelector.Items.Add("tap_down");
            gestureSelector.Items.Add("tap_front");
            gestureSelector.Items.Add("tap_back");
            gestureSelector.Items.Add("scratchout");
            gestureSelector.Items.Add("horizontal_square");
            gestureSelector.Items.Add("vertical_square");
            gestureSelector.Items.Add("horizontal_circle");
            gestureSelector.Items.Add("vertical_circle");         
            gestureSelector.Items.Add("horizontal_x");
            gestureSelector.Items.Add("vertical_x");
            gestureSelector.Items.Add("horizontal_L");
            gestureSelector.Items.Add("vertical_L");
            //gestureSelector.Items.Add("z");
            gestureSelector.Items.Add("detect");            
        }

        private void use3DGestures_Unchecked(object sender, RoutedEventArgs e)
        {
            gestureSelector.Items.Clear();
            gestureSelector.Items.Add("swipe_left");
            gestureSelector.Items.Add("swipe_right");
            gestureSelector.Items.Add("swipe_up");
            gestureSelector.Items.Add("swipe_down");
            gestureSelector.Items.Add("tap_left");
            gestureSelector.Items.Add("tap_right");
            gestureSelector.Items.Add("tap_up");
            gestureSelector.Items.Add("tap_down");
            gestureSelector.Items.Add("scratchout");
            gestureSelector.Items.Add("circle");
            gestureSelector.Items.Add("square");
            gestureSelector.Items.Add("x");
            gestureSelector.Items.Add("L");
            //gestureSelector.Items.Add("c");
            //gestureSelector.Items.Add("two_handed_fb");
            //gestureSelector.Items.Add("two_handed_lr");
            gestureSelector.Items.Add("detect");
        }
        
    }
}
