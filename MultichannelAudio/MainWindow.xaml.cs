using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using NAudio;
using NAudio.Wave;
using NAudio.Dsp;
using Exocortex.DSP;
using Microsoft.Research.DynamicDataDisplay;
using Microsoft.Research.DynamicDataDisplay.DataSources;

namespace MultichannelAudio
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        //private WaveOut waveOut;
        private WasapiOut wOut;
        private WaveIn waveIn;

        public int waveOutChannels;        

        public int buffersize = 2048;
        public int[] bin;
        public float[] sampledata;
        public float[] inbetween;
        bool init_inbetween = true;
        ComplexF[] indata;

        double[] filteredindata;
        double[] priori;

        EnumerableDataSource<int> bins;
        EnumerableDataSource<float> rawIn;
        EnumerableDataSource<double> fftIn;

        int[] channelLabel;
        int[] velocity;
        Rectangle[] bars;
        
        //EnumerableDataSource<int> channelSource;
        //EnumerableDataSource<int> velocitySource;

        int selectedChannels = 1;
        List<int> frequencies;
        List<int> centerbins;

        List<KeyFrequency> KF;

        StreamWriter svmOutStream = null;
        //List<StreamWriter> Foutstream;

        private bool graphing = true;
        private bool recording = false;        

        public MainWindow()
        {         
            InitializeComponent();
            this.WindowStartupLocation = System.Windows.WindowStartupLocation.CenterScreen;
            int waveInDevices = WaveIn.DeviceCount;            
            for (int waveInDevice = 0; waveInDevice < waveInDevices; waveInDevice++)
            {
                WaveInCapabilities deviceInfo = WaveIn.GetCapabilities(waveInDevice);                
                Console.WriteLine("Device {0}: {1}, {2} channels",
                    waveInDevice, deviceInfo.ProductName, deviceInfo.Channels);           
            }

            this.gestureSelector.SelectedIndex = 0;

            waveIn = new WaveIn();
            waveIn.BufferMilliseconds = 47*buffersize/2048;
            waveIn.DeviceNumber = 0;
            waveIn.WaveFormat = new WaveFormat(44100, 32, 1);
            waveIn.DataAvailable += waveIn_DataAvailable;

            try
            {
                waveIn.StartRecording();
            }
            catch(NAudio.MmException e)
            {
                Console.WriteLine(e.ToString() + "\nPlug in a microphone!");
            }

            bin = new int[buffersize * 2];
            sampledata = new float[buffersize * 2];
            priori = new double[buffersize * 2];
            
            channelLabel = new int[1];
            channelLabel[0] = 1;
            velocity = new int[1];
            velocity[0] = 0;
            bars = new Rectangle[1];
            bars[0] = createBar(channelLabel[0], selectedChannels, velocity[0], 33);
            for (int i = 0; i < buffersize * 2; i++)
            {
                bin[i] = i;
                sampledata[i] = 0;
                priori[i] = 0;
                
            }

            chart1.Viewport.Visible = new DataRect(0, -1.0, buffersize * 2, 2.0);
            
            //chart2.Viewport.Visible = new DataRect(1620, -50, 280,110);
            chart2.Viewport.Visible = new DataRect(1620, 0, 280, 110);   
                     
            //chart3.Viewport.Visible = new DataRect(1, -30, 1, 60);
            
            bins = new EnumerableDataSource<int>(bin);
            bins.SetXMapping(x => x);

            rawIn = new EnumerableDataSource<float>(sampledata);
            rawIn.SetYMapping(y => y);

            CompositeDataSource comp1 = new CompositeDataSource(bins, rawIn);
            chart1.AddLineGraph(comp1);

            CompositeDataSource comp2 = new CompositeDataSource(bins, rawIn);
            chart2.AddLineGraph(comp2);

            //channelSource = new EnumerableDataSource<int>(channelLabel);
            //channelSource.SetXMapping(x => x);
            //velocitySource = new EnumerableDataSource<int>(velocity);
            //velocitySource.SetYMapping(y => y);
            //CompositeDataSource comp3 = new CompositeDataSource(channelSource, velocitySource);
            //chart3.AddLineGraph(comp3);
            _barcanvas.Background = new SolidColorBrush(Colors.LightBlue);
        }

        void waveIn_DataAvailable(object sender, WaveInEventArgs e)
        {
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
                                        (e.Buffer[index*4 + 0]));
                
                float sample32 = sample / 2147483648f;
                
                if (index >= (buffersize - inbetween.Length))
                    sampledata[index] = inbetween[index - buffersize + inbetween.Length];
                else
                    sampledata[index] = sampledata[index+buffersize+inbetween.Length]; 
                sampledata[index+buffersize] = sample32;
               
            }

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
            if (graphing)
            {
                bins = new EnumerableDataSource<int>(bin);
                bins.SetXMapping(x => x);

                rawIn = new EnumerableDataSource<float>(sampledata);
                rawIn.SetYMapping(y => y);
                chart1.Children.RemoveAll<LineGraph>();
                CompositeDataSource comp1 = new CompositeDataSource(bins, rawIn);
                chart1.AddLineGraph(comp1, Colors.Red);

                //fftIn = new EnumerableDataSource<Complex>(indata);
                //fftIn.SetYMapping(y => Math.Sqrt(y.X*y.X + y.Y*y.Y));                
                //fftIn.SetYMapping(y => Math.Sqrt(y.Re * y.Re + y.Im * y.Im));

                fftIn = new EnumerableDataSource<double>(filteredindata);
                fftIn.SetYMapping(y => y);

                //fftIn = new EnumerableDataSource<ComplexF>(indata);
                //fftIn.SetYMapping(y => mag2db(y));


                chart2.Children.RemoveAll<LineGraph>();
                CompositeDataSource comp2 = new CompositeDataSource(bins, fftIn);
                chart2.AddLineGraph(comp2, Colors.Red);
            }

            if ((wOut != null) && recording)
            {
                //for (int i = 0; i < frequencies.Count; i++)
                //{
                //    if (Foutstream.Count < i+1)
                //    {
                //        Foutstream.Add(new StreamWriter(filenameBox.Text + i.ToString() + ".dat"));
                //        Console.WriteLine(filenameBox.Text + i.ToString() + ".dat");
                //    }
                //}

                //channelSource = new EnumerableDataSource<int>(channelLabel);

                _barcanvas.Children.Clear();

                if(svmOutStream == null)
                    svmOutStream = new StreamWriter("test.ds");

                KF = new List<KeyFrequency>();
                gestureDetected.Text = "";
                for (int i = 0; i < frequencies.Count; i++)
                {                    
                    KF.Add(new KeyFrequency(frequencies.ElementAt(i), i + 1, 33, filteredindata, centerbins.ElementAt(i)));
                    //Foutstream.ElementAt(i).WriteLine(KF.ElementAt(i));
                    svmOutStream.WriteLine(gestureSelector.SelectedIndex + " " + KF.ElementAt(i));
                    velocity[i] = KF.ElementAt(i).state;
                    //Console.WriteLine(velocity[i]);
                    gestureDetected.Text += " " + KF.ElementAt(i).state.ToString();
                    bars[i] = createBar(channelLabel[i], selectedChannels, velocity[i], 33);
                    _barcanvas.Children.Add(bars[i]);
                }

                //chart3.Children.RemoveAll<LineGraph>();

                //channelSource = new EnumerableDataSource<int>(channelLabel);
                //channelSource.SetXMapping(x => x);
                //velocitySource = new EnumerableDataSource<int>(velocity);
                //velocitySource.SetYMapping(y => y);
                //CompositeDataSource comp3 = new CompositeDataSource(channelSource, velocitySource);
                //chart3.AddLineGraph(comp3, Colors.Red);       


            }


        }

        public void bufferFFT()
        {
            indata = new ComplexF[buffersize * 2];
            for (int i = 0; i < buffersize * 2; i++)
            {
                indata[i].Re = sampledata[i] *(float)FastFourierTransform.HammingWindow(i, buffersize * 2);
                indata[i].Im = 0;
            }
            //FastFourierTransform.FFT(true, 11, indata);
            Exocortex.DSP.Fourier.FFT(indata, buffersize * 2, Exocortex.DSP.FourierDirection.Forward);
            filteredindata = filterMean(indata,1.3);
        }

        private void button1_Click(object sender, EventArgs e)
        {
            StartStopSineWave();
        }

        private void record_Click(object sender, EventArgs e)
        {
            if (recording)
            {
                recording = false;
                //for (int i = 0; i < Foutstream.Count; i++)
                //{                   
                //    Foutstream.ElementAt(i).Close();                    
                //}                
                button2.Content = "Record";
            }
            else
            {
                recording = true;
                button2.Content = "Stop";
            }
        }


        private void StartStopSineWave()
        {
            if (wOut == null)
            {
                button1.Content = "Stop Sound";
                //string str = channelSelector.Text;
                //selectedChannels = Int32.Parse(str);
                Console.WriteLine("User Selected Channels: " + selectedChannels);
                WaveOutCapabilities outdeviceInfo = WaveOut.GetCapabilities(0);                
                waveOutChannels = outdeviceInfo.Channels;
                wOut = new WasapiOut(NAudio.CoreAudioApi.AudioClientShareMode.Shared, 0);
                int waveOutDevices = WaveOut.DeviceCount;
                for(int i = 0; i< waveOutDevices; i++)
                {
                    outdeviceInfo = WaveOut.GetCapabilities(i);
                    Console.WriteLine("Device {0}: {1}, {2} channels",
                            i, outdeviceInfo.ProductName, outdeviceInfo.Channels);
                }

                /*
                var sineWaveProvider = new SineWaveProvider32();
                sineWaveProvider.SetWaveFormat(44100, 1); // 44.1kHz mono
                sineWaveProvider.Frequency = 19000;
                sineWaveProvider.Amplitude = 0.25f;
                var sineWaveProvider2 = new SineWaveProvider32();
                sineWaveProvider2.SetWaveFormat(44100, 1); // 44.1kHz mono
                sineWaveProvider2.Frequency = 18000;
                sineWaveProvider2.Amplitude = 0.25f;
                
                List<IWaveProvider> inputs = new List<IWaveProvider>();
                inputs.Add(sineWaveProvider);
                inputs.Add(sineWaveProvider2);
                */

                List<IWaveProvider> inputs = new List<IWaveProvider>();
                frequencies = new List<int>();
                centerbins = new List<int>();
                //Foutstream = new List<StreamWriter>();
                for (int c = 0; c < selectedChannels; c++)
                {
                    inputs.Add(new SineWaveProvider32(18000+c*700, 0.25f, 44100, 1));
                    frequencies.Add(18000+c*700);
                    centerbins.Add((int)Math.Round((18000 + c * 700) / 10.768));                  
                    
                    Console.WriteLine(centerbins.ElementAt(c));                    
                }

                var splitter = new MultiplexingWaveProvider(inputs, selectedChannels);
                try
                {
                    wOut.Init(splitter);
                    wOut.Play();
                }
                catch(System.ArgumentException)
                {
                    Console.WriteLine("Invalid audio channel count. Please select a lower number of audio channels");
                }
                //waveOut = new WaveOut();
                //waveOut.Init(sineWaveProvider);                    
                //waveOut.Init(splitter);
                
                Console.WriteLine("Number of Channels: " + wOut.OutputWaveFormat.Channels);
            }
            else
            {
                wOut.Stop();
                wOut.Dispose();
                wOut = null;
                button1.Content = "Generate Sound";

                frequencies.Clear();
                centerbins.Clear();
                //Foutstream.Clear();
            }
        }

        public float mag2db(ComplexF y)
        {
            return 20.0f * (float)Math.Log10(Math.Sqrt(y.Re * y.Re + y.Im * y.Im)/.02);
        }

        public double[] filterMean(ComplexF[] data, double factor)
        {
            double[] outdata = new double[data.Length];
            double min = Double.PositiveInfinity;
            double mean = 0;
            for (int i = 0; i < data.Length; i++)
            {
                outdata[i] = mag2db(data[i]);
                min = Math.Min(outdata[i], min);
            }
            /*
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
            }*/
            //Console.WriteLine("Mean: " + mean + " Min: " + min);
            return outdata;
        }

        public double[] filterMovingAverage(ComplexF[] input, int width)
        {
            double[] outdata = new double[input.Length];
            if (width == 1)
                for (int i = 0; i < input.Length; i++)
                    outdata[i] = mag2db(input[i]);
                
            if (width == 3)
                for (int i = 1; i < input.Length - 1; i++)            
                    outdata[i] = (mag2db(input[i - 1]) + mag2db(input[i]) + mag2db(input[i + 1])) / 3;
            
            if (width == 5)
                for (int i = 2; i < input.Length - 2; i++)            
                    outdata[i] = (mag2db(input[i - 2]) + mag2db(input[i - 1]) + mag2db(input[i]) + mag2db(input[i + 1]) + mag2db(input[i + 2])) / 5;
            
            return outdata;
        }

        public double[] filterHighPass(ComplexF[] input)
        {
            double[] outdata = new double[input.Length];
            outdata[0] = mag2db(input[0]);
            for (int i = 1; i < input.Length; i++)
            {
                outdata[i] = mag2db(input[i]) - mag2db(input[i - 1]);
            }
            return outdata;
        }

        public double[] filterGaussian(ComplexF[] input, int sigma)
        {
            double[] outdata = new double[input.Length];
            double[] converted = new double[input.Length];
            for (int i = 0; i < input.Length; i++)
            {
                converted[i] = mag2db(input[i]);
                outdata[i] = mag2db(input[i]);
            }

            double constant = 1 / (Math.Sqrt(2 * Math.PI) * sigma);
            double[] window = new double[1 + 2 * sigma];
            for (int i = 0; i < window.Length; i++)
            {
                window[i] = constant * Math.Exp(-1 * Math.Pow(i - sigma, 2) / 2 * Math.Pow(sigma, 2));
            }
            for (int i = sigma; i < input.Length - sigma; i++)
            {
                outdata[i] = 0;
                for (int j = -sigma; j <= sigma; j++)
                {
                    outdata[i] += outdata[i + j] * window[j + sigma];
                }
            }
            return outdata;
        }

        protected override void OnClosed(EventArgs e)
        {
            base.OnClosed(e);
            if (svmOutStream != null)
            {
                svmOutStream.Close();
                svmOutStream.Dispose();
            }
        }
        
        private void channelSelector_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            selectedChannels = (sender as ComboBox).SelectedIndex + 1;
            //selectedChannels = Int32.Parse(str);
            channelLabel = new int[selectedChannels];
            velocity = new int[selectedChannels];
            bars = new Rectangle[selectedChannels];
            //chart3.Viewport.Visible = new DataRect(1, -30, selectedChannels, 60);

            //channelSource = new EnumerableDataSource<int>(channelLabel);
            //channelSource.SetXMapping(x => x);

            for (int i = 0; i < selectedChannels; i++)
            {
                channelLabel[i] = i + 1;
                velocity[i] = 0;
                bars[i] = createBar(i+1, selectedChannels, 0, 33);
            }

        }

        private Rectangle createBar(int channel, int cmax, int velocity, int vmax)
        {
            Rectangle outRect = new Rectangle();
            outRect.Width = _barcanvas.ActualWidth / cmax;
            outRect.Height = (_barcanvas.ActualHeight*Math.Abs(velocity))/(2*vmax);

            Canvas.SetLeft(outRect,(outRect.Width*(channel-1)));

            if (velocity < 0)
            {
                Canvas.SetTop(outRect, _barcanvas.ActualHeight / 2);
                outRect.Stroke = new SolidColorBrush(Colors.Black);
                outRect.Fill = new SolidColorBrush(Colors.Red);
            }
            else
            {
                Canvas.SetTop(outRect, (_barcanvas.ActualHeight / 2)-((_barcanvas.ActualHeight * Math.Abs(velocity)) / (2 * vmax)));
                outRect.Stroke = new SolidColorBrush(Colors.Black);
                outRect.Fill = new SolidColorBrush(Colors.Green);
            }

            Console.WriteLine(outRect.Width + " " + outRect.Height);
            return outRect;
        }

    }
}
