using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
using GestureTests.Gesture;
using GestureTests.Util;
using GestureTests.Types;

/*
 
Author: Salman Cheema
University of Central Florida
 
Email: salmanc@cs.ucf.edu
 
Released as part of the 3D Gesture Database analysed in
 
"Salman Cheema, Michael Hoffman, Joseph J. LaViola Jr., 3D Gesture classification with linear acceleration and angular velocity 
sensing devices for video games, Entertainment Computing, Volume 4, Issue 1, February 2013, Pages 11-24, ISSN 1875-9521, 10.1016/j.entcom.2012.09.002"
 
*/

namespace GestureTests.Data
{
    /// <summary>
    /// Dataset for a single user. Contains 3 different subsets:
    /// Training samples (625), Correctly Classified Gameplay samples (varies for each user), Incorrectly Classified Gameplay samples (varies for each user)
    /// </summary>
    public class UserDataSet
    {
        /// <summary>
        /// Path where this dataset is stored on disk.
        /// </summary>
        public string Path;

        /// <summary>
        /// Set of training data. 25 samples for each gesture = 625 training samples/user.
        /// </summary>
        public List<GestureSample> TrainingSamples { get; protected set; }

        /// <summary>
        /// Set of gestures that were correctly classified during gameplay.
        /// </summary>
        public List<GestureSample> CorrectlyClassifiedGameplaySamples { get; protected set; }

        /// <summary>
        /// Set of gestures that were incorrectly classified during gameplay.
        /// </summary>
        public List<GestureSample> IncorrectlyClassifiedGameplaySamples { get; protected set; }

        public UserDataSet(string path)
        {
            try
            {
                ///load the dataset for this user.
                this.Path = path;
                this.TrainingSamples = new List<GestureSample>();
                string[] training_files = Directory.GetFiles(this.Path);

                foreach (string tFile in training_files)
                {
                    GestureSample sample = LoadSample(tFile);
                    this.TrainingSamples.Add(sample);
                }
                /*
                //load training data
                this.TrainingSamples = new List<GestureSample>();
                string training_path = path + "\\training\\";
                string[] training_files = Directory.GetFiles(training_path);

                foreach (string tFile in training_files)
                {
                    GestureSample sample = LoadSample(tFile);
                    this.TrainingSamples.Add(sample);
                }
                */
                //go through gameplay runs and load correct/incorrect gameplay samples
                this.IncorrectlyClassifiedGameplaySamples = new List<GestureSample>();
                this.CorrectlyClassifiedGameplaySamples = new List<GestureSample>();
                /*
                string gameplay_path = path + "\\gameplay";
                int count = 0;
                while (Directory.Exists(gameplay_path + count))
                {
                    string ith_gameplay_path = gameplay_path + count + "\\";
                    string[] gameplay_samples = Directory.GetFiles(ith_gameplay_path);
                    foreach (string gFile in gameplay_samples)
                    {
                        if (gFile.Contains(".txt")) continue;

                        GestureSample sample = LoadSample(gFile);

                        if (gFile.Contains("-incorrect"))
                            this.IncorrectlyClassifiedGameplaySamples.Add(sample);
                        else
                            this.CorrectlyClassifiedGameplaySamples.Add(sample);
                    }
                    ++count;
                }
                */
                Console.WriteLine("Found " + this.TrainingSamples.Count + " Training, " + this.CorrectlyClassifiedGameplaySamples.Count + " Correct Gameplay, " + this.IncorrectlyClassifiedGameplaySamples.Count + " Incorrect Gameplay at [" + path + "]");
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                Console.WriteLine(e.StackTrace);
            }
        }

        /// <summary>
        /// Parses a given file and loads the gesture data from it.
        /// </summary>
        /// <param name="filename">path of file containing an instance of a gesture</param>
        /// <returns>gesture sample loaded from the file. Returns null if an I/O error occurs.</returns>
        private GestureSample LoadSample(string filename)
        {
            if (Config.Use3DMode == true)
                return Load3DSample(filename);
            GestureSample sample = null;

            GestureType gesture = GestureType.unknown;
            float duration = float.NaN;
            bool rightHanded = false;
            List<Vector2> interpretedPoints = new List<Vector2>();
            List<Vector2> velocities = new List<Vector2>();
            List<Vector2> inverseVelocities = new List<Vector2>();
            List<Vector2> strokePoints = new List<Vector2>();
            List<float> angles = new List<float>();            
            try
            {
                StreamReader reader = File.OpenText(filename);

                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();

                    //end of file
                    if (line == null) continue;

                    //skip lines starting with '#' are comments
                    if (line.StartsWith("#")) continue;

                    //ignore empty lines
                    if (line == "") continue;

                    //split the line by the 'space' character
                    string[] tokens = line.Split(" ".ToCharArray());

                    //the first token provides information about the type of data to follow
                    switch (tokens[0])
                    {
                        case "GestureName:":
                            gesture = ReadGestureType(tokens[1]);
                            break;
                        case "Duration(ms):":
                            duration = float.Parse(tokens[1]);
                            break;
                        case "Handedness:":
                            rightHanded = (tokens[1] == "right");
                            break;
                        case "SpeakerAngles:":
                            int numAngles = int.Parse(tokens[1]);
                            for (int i = 0; i < numAngles; ++i)
                            {
                                string angle = reader.ReadLine();
                                string[] theta = angle.Split(" ,".ToCharArray());
                                angles.Add(float.Parse(theta[0]));
                            }

                            break;                        
                        case "InterpretedPoints:":
                            int numPoints = int.Parse(tokens[1]);

                            //read the points from succeeding lines.
                            for (int i = 0; i < numPoints; ++i)
                            {
                                string point = reader.ReadLine();
                                string[] xy = point.Split(" ,".ToCharArray());
                                interpretedPoints.Add(new Vector2(float.Parse(xy[0]), float.Parse(xy[1])));
                            }

                            break;
                        case "StrokePoints:":
                            int numStrokePoints = int.Parse(tokens[1]);

                            //read datapoints from succeeding lines.
                            for (int i = 0; i < numStrokePoints; ++i)
                            {
                                string point = reader.ReadLine();
                                string[] xy = point.Split(" ,".ToCharArray());
                                //180f, 134.34f
                                strokePoints.Add(new Vector2(float.Parse(xy[0]), float.Parse(xy[1])));
                            }
                            break;

                        case "Velocities:":
                            int numVelocities = int.Parse(tokens[1]);

                            //read datapoints from succeeding lines.
                            for (int i = 0; i < numVelocities; ++i)
                            {
                                string point = reader.ReadLine();
                                string[] xy = point.Split(" ,".ToCharArray());
                                velocities.Add(new Vector2(float.Parse(xy[0]), float.Parse(xy[1])));
                            }
                            break;
                        case "InverseVelocities:":
                            int numInverseVelocities = int.Parse(tokens[1]);

                            //read datapoints from succeeding lines.
                            for (int i = 0; i < numInverseVelocities; ++i)
                            {
                                string point = reader.ReadLine();
                                string[] xy = point.Split(" ,".ToCharArray());
                                inverseVelocities.Add(new Vector2(float.Parse(xy[0]), float.Parse(xy[1])));
                            }
                            break;
                    }
                }
                if (strokePoints.Count == 0)
                    strokePoints = GenerateStroke(interpretedPoints);
                reader.Close();
                sample = new GestureSample(gesture, rightHanded, duration, angles, interpretedPoints, strokePoints, velocities, inverseVelocities);
                sample.ComputeFeatures(Config.FeaturesToUse);
            }
            catch (Exception e)
            {

                Console.WriteLine(e.Message);
                Console.WriteLine(e.StackTrace);
            }


            return sample;
        }

        private GestureSample Load3DSample(string filename)
        {
            GestureSample sample = null;

            GestureType gesture = GestureType.unknown;
            float duration = float.NaN;
            bool rightHanded = false;
            List<Vector3> interpretedPoints = new List<Vector3>();
            List<Vector3> velocities = new List<Vector3>();
            List<Vector3> inverseVelocities = new List<Vector3>();
            List<Vector3> strokePoints = new List<Vector3>();
            List<float> angles = new List<float>();
            List<float> elevations = new List<float>();

            try
            {
                StreamReader reader = File.OpenText(filename);

                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();

                    //end of file
                    if (line == null) continue;

                    //skip lines starting with '#' are comments
                    if (line.StartsWith("#")) continue;

                    //ignore empty lines
                    if (line == "") continue;

                    //split the line by the 'space' character
                    string[] tokens = line.Split(" ".ToCharArray());

                    //the first token provides information about the type of data to follow
                    switch (tokens[0])
                    {
                        case "GestureName:":
                            gesture = ReadGestureType(tokens[1]);
                            break;
                        case "Duration(ms):":
                            duration = float.Parse(tokens[1]);
                            break;
                        case "Handedness:":
                            rightHanded = (tokens[1] == "right");
                            break;
                        case "SpeakerAngles:":
                            int numAngles = int.Parse(tokens[1]);
                            for (int i = 0; i < numAngles; ++i)
                            {
                                string angle = reader.ReadLine();
                                string[] theta = angle.Split(" ,".ToCharArray());
                                angles.Add(float.Parse(theta[0]));
                            }

                            break;
                        case "SpeakerElevations:":
                            int numEle = int.Parse(tokens[1]);
                            for (int i = 0; i < numEle; ++i)
                            {
                                string angle = reader.ReadLine();
                                string[] ele = angle.Split(" ,".ToCharArray());
                                elevations.Add(float.Parse(ele[0]));
                            }

                            break;
                        case "InterpretedPoints:":
                            int numPoints = int.Parse(tokens[1]);

                            //read the points from succeeding lines.
                            for (int i = 0; i < numPoints; ++i)
                            {
                                string point = reader.ReadLine();
                                string[] xyz = point.Split(" ,".ToCharArray());
                                interpretedPoints.Add(new Vector3(float.Parse(xyz[0]), float.Parse(xyz[1]), float.Parse(xyz[2])));
                            }

                            break;
                        case "StrokePoints:":
                            int numStrokePoints = int.Parse(tokens[1]);

                            //read datapoints from succeeding lines.
                            for (int i = 0; i < numStrokePoints; ++i)
                            {
                                string point = reader.ReadLine();
                                string[] xyz = point.Split(" ,".ToCharArray());
                                strokePoints.Add(new Vector3(float.Parse(xyz[0]), float.Parse(xyz[1]), float.Parse(xyz[2])));
                            }
                            break;

                        case "Velocities:":
                            int numVelocities = int.Parse(tokens[1]);

                            //read datapoints from succeeding lines.
                            for (int i = 0; i < numVelocities; ++i)
                            {
                                string point = reader.ReadLine();
                                string[] xyz = point.Split(" ,".ToCharArray());
                                velocities.Add(new Vector3(float.Parse(xyz[0]), float.Parse(xyz[1]), float.Parse(xyz[2])));
                            }
                            break;
                        case "InverseVelocities:":
                            int numInverseVelocities = int.Parse(tokens[1]);

                            //read datapoints from succeeding lines.
                            for (int i = 0; i < numInverseVelocities; ++i)
                            {
                                string point = reader.ReadLine();
                                string[] xyz = point.Split(" ,".ToCharArray());
                                inverseVelocities.Add(new Vector3(float.Parse(xyz[0]), float.Parse(xyz[1]), float.Parse(xyz[2])));
                            }
                            break;
                    }
                }
                if (strokePoints.Count == 0)
                    strokePoints = Generate3DStroke(interpretedPoints);
                reader.Close();
                sample = new GestureSample(gesture, rightHanded, duration, angles, elevations, interpretedPoints, strokePoints, velocities, inverseVelocities);
                sample.ComputeFeatures(Config.FeaturesToUse);
            }
            catch (Exception e)
            {

                Console.WriteLine(e.Message);
                Console.WriteLine(e.StackTrace);
            }


            return sample;
        }

        private List<Vector2> GenerateStroke(List<Vector2> interpretedPoints)
        {
            List<Vector2> S = new List<Vector2>();
            S.Add(new Vector2(180f, 134.34f));
            for (int i = 0; i < interpretedPoints.Count; i++)
            {
                S.Add(new Vector2(S[i].X - interpretedPoints[i].X, S[i].Y - interpretedPoints[i].Y));
            }

            return S;
        }

        private List<Vector3> Generate3DStroke(List<Vector3> interpretedPoints)
        {
            List<Vector3> S = new List<Vector3>();
            S.Add(new Vector3(0,0,0));
            for (int i = 0; i < interpretedPoints.Count; i++)
            {
                S.Add(new Vector3(S[i].X + interpretedPoints[i].X, S[i].Y + interpretedPoints[i].Y, S[i].Z + interpretedPoints[i].Z));
            }

            return S;
        }

        /// <summary>
        /// Converts a gesture name in text to a 'GestureType'.
        /// </summary>
        /// <param name="gestureName"></param>
        /// <returns></returns>
        private GestureType ReadGestureType(string gestureName)
        {
            switch (gestureName)
            {
                case "swipe_left": return GestureType.swipe_left;
                case "swipe_right": return GestureType.swipe_right;
                case "swipe_up": return GestureType.swipe_up;
                case "swipe_down": return GestureType.swipe_down;
                case "swipe_front": return GestureType.swipe_front;
                case "swipe_back": return GestureType.swipe_back;
                case "tap_left": return GestureType.tap_left;
                case "tap_right": return GestureType.tap_right;
                case "tap_up": return GestureType.tap_up;
                case "tap_down": return GestureType.tap_down;
                case "tap_front": return GestureType.tap_front;
                case "tap_back": return GestureType.tap_back;
                case "scratchout": return GestureType.scratchout;
                case "square": return GestureType.square;
                case "x": return GestureType.x;
                case "c": return GestureType.c;
                case "circle": return GestureType.circle;
                case "two_handed_fb": return GestureType.two_handed_fb;
                case "two_handed_lr": return GestureType.two_handed_lr;
                case "horizontal_circle": return GestureType.horizontal_circle;
                case "vertical_circle": return GestureType.vertical_circle;
                case "spiral": return GestureType.spiral;
                case "arm_lift": return GestureType.arm_lift;                
                case "triangle": return GestureType.triangle;
                case "z": return GestureType.z;                                
                case "side_x": return GestureType.side_x;
                case "rectangle": return GestureType.rectangle;                
                case "check": return GestureType.check;
                case "caret": return GestureType.caret;
                case "zigzag": return GestureType.zigzag;
                case "arrow": return GestureType.arrow;
                case "left_bracket": return GestureType.left_bracket;
                case "right_bracket": return GestureType.right_bracket;
                case "star": return GestureType.star;                
            }

            return GestureType.unknown;
        }     
    }
}
