using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using GestureTests.Gesture;

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
    /// Class to load the dataset. This class can also export the given dataset as '.arff' files for use in WEKA.
    /// </summary>
    public class DataLoader
    {
        /// <summary>
        /// Loads the gesture dataset from a given location.
        /// </summary>
        /// <param name="dataPath">path from which to load data.</param>
        /// <returns>returns the dataset. returns empty list if no data found or some error occurs.</returns>
        public static List<UserDataSet> LoadGestureDataFrom(string dataPath)
        {
            List<UserDataSet> gestureData = new List<UserDataSet>();
            try
            {
                if (!Directory.Exists(dataPath))
                {
                    Console.WriteLine("Unable to find User Data Directory");
                    return gestureData;
                }

                string[] usernames = Directory.GetDirectories(dataPath);
                if (usernames.Length == 0)
                {
                    Console.WriteLine("Found NO User Data at [" + dataPath + "]");
                    return gestureData;
                }

                Console.WriteLine("Found Data for " + usernames.Length + " participants");

                for (int i = 0; i < usernames.Length; ++i)
                    gestureData.Add(new UserDataSet(usernames[i]));
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
            return gestureData;
        }

        /// <summary>
        /// Loads a given dataset, and exports each gesture sample's feature vector as '.arff' files to be analyzed by WEKA.
        /// For each user, 3 '.arff' files are generated, pertaining to training data, correctly classifier gameplay data, and incorrently classified gameplay data.
        /// </summary>
        /// <param name="dataPath">Source path for dataset.</param>
        /// <param name="outputPath">Target path for exporting dataset.</param>
        public static void ExportFeaturesAsARFF(string dataPath, string outputPath)
        {

            //load the experiment data
            List<UserDataSet> data = LoadGestureDataFrom(dataPath);

            foreach (UserDataSet data_i in data)
            {
                int index = data_i.Path.LastIndexOf("\\");
                string name = data_i.Path.Substring(index + 1, data_i.Path.Length - index - 1);

                //export the three different lists of gestures for each user in a separate file.
                string filename_prefix = outputPath + name + "-";
                WriteARFF(name, filename_prefix + "training.arff", data_i.TrainingSamples);
                //WriteARFF(name, filename_prefix + "correct_gameplay.arff", data_i.CorrectlyClassifiedGameplaySamples);
                //WriteARFF(name, filename_prefix + "incorrect_gameplay.arff", data_i.IncorrectlyClassifiedGameplaySamples);
            }

        }

        /// <summary>
        /// Writes the feature vectors for a list of gestures to an '.arff' file.
        /// </summary>
        /// <param name="userid">username for user who recorded the gestures.</param>
        /// <param name="filename">path of output file.</param>
        /// <param name="samples">list of gestures.</param>
        private static void WriteARFF(string userid, string filename, List<GestureSample> samples)
        {
            if (Config.Use3DMode)
            {
                WriteARFF3D(userid, filename, samples);
                return;
            }
            //create file
            StreamWriter file = File.CreateText(filename);

            //write arff header
            file.Write("@RELATION " + filename + "\r\r");

            //file.Write("@ATTRIBUTE username STRING\r");
            file.Write("@ATTRIBUTE GestureType  {swipe_left,swipe_right,swipe_up,swipe_down,tap_left,tap_right,tap_up,tap_down,scratchout,circle,square,x,c,two_handed_fb,two_handed_lr,unknown}\r");
            //file.Write("@ATTRIBUTE Angle1 NUMERIC\r");
            //file.Write("@ATTRIBUTE Angle2 NUMERIC\r");
            file.Write("@ATTRIBUTE MinX NUMERIC\r");
            file.Write("@ATTRIBUTE MinY NUMERIC\r");
            file.Write("@ATTRIBUTE MaxX NUMERIC\r");
            file.Write("@ATTRIBUTE MaxY NUMERIC\r");            
            file.Write("@ATTRIBUTE MeanX NUMERIC\r");
            file.Write("@ATTRIBUTE MeanY NUMERIC\r");            
            file.Write("@ATTRIBUTE MedianX NUMERIC\r");
            file.Write("@ATTRIBUTE MedianY NUMERIC\r");
            file.Write("@ATTRIBUTE DiagonalLength NUMERIC\r");
            file.Write("@ATTRIBUTE StartAngleSinXY NUMERIC\r");
            file.Write("@ATTRIBUTE StartAngleCosXY NUMERIC\r");
            file.Write("@ATTRIBUTE FirstLastAngleSinXY NUMERIC\r");
            file.Write("@ATTRIBUTE FirstLastAngleCosXY NUMERIC\r");
            file.Write("@ATTRIBUTE TotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE AbsoluteTotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE SquaredTotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE FirstLastDistance NUMERIC\r");
            file.Write("@ATTRIBUTE TotalDistance NUMERIC\r");
            file.Write("@ATTRIBUTE MaxAccelerationSquared NUMERIC\r");
            file.Write("@ATTRIBUTE ZeroCrossingsX NUMERIC\r");
            file.Write("@ATTRIBUTE ZeroCrossingsY NUMERIC\r");
            /*file.Write("@ATTRIBUTE Inverse_MinX NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_MinY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_MaxX NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_MaxY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_MeanX NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_MeanY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_MedianX NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_MedianY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_DiagonalLength NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_StartAngleSinXY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_StartAngleCosXY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_FirstLastAngleSinXY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_FirstLastAngleCosXY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_TotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_AbsoluteTotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_SquaredTotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_FirstLastDistance NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_TotalDistance NUMERIC\r");
            file.Write("@ATTRIBUTE Inverse_MaxAccelerationSquared NUMERIC\r");*/
            /*file.Write("@ATTRIBUTE Velocity_MinX NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MinY NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MaxX NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MaxY NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MeanX NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MeanY NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MedianX NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MedianY NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MinX NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MinY NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MaxX NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MaxY NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MeanX NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MeanY NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MedianX NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MedianY NUMERIC\r");*/         
            file.Write("@ATTRIBUTE Stroke_MinX NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MinY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MaxX NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MaxY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MeanX NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MeanY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MedianX NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MedianY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_DiagonalLength NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_StartAngleSinXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_StartAngleCosXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_FirstLastAngleSinXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_FirstLastAngleCosXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_TotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_AbsoluteTotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_SquaredTotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_FirstLastDistance NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_TotalDistance NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MaxAccelerationSquared NUMERIC\r");
            /*file.Write("@ATTRIBUTE Stroke_ZeroCrossingsX NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_ZeroCrossingsY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_BoundingBoxHeight NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_BoundingBoxWidth NUMERIC\r");*/
            file.Write("@ATTRIBUTE Stroke_SideRatioFirst NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_SideRatioLast NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_TopRatioFirst NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_TopRatioLast NUMERIC\r");
            /*file.Write("@ATTRIBUTE Stroke_AspectRatio NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_Regression NUMERIC\r");*/
            file.Write("@ATTRIBUTE Duration NUMERIC\r");


            //write out the feature vector for this sample
            file.Write("\r@DATA\r");
            foreach (GestureSample gs in samples)
            {
                string datum = gs.Gesture.ToString();
                for (int i = 0; i < gs.Features.Length; ++i)
                    datum += "," + gs.Features[i];
                file.Write(datum + "\r");
            }

            file.Flush();
            file.Close();
        }

        private static void WriteARFF3D(string userid, string filename, List<GestureSample> samples)
        {
            //create file
            StreamWriter file = File.CreateText(filename);

            //write arff header
            file.Write("@RELATION " + filename + "\r\r");

            //file.Write("@ATTRIBUTE username STRING\r");
            file.Write("@ATTRIBUTE GestureType  {swipe_left,swipe_right,swipe_up,swipe_down,swipe_front,swipe_back,tap_left,tap_right,tap_up,tap_down,tap_front,tap_back,scratchout,circle,square,x,c,two_handed_fb,two_handed_lr,horizontal_circle,vertical_circle,spiral,arm_lift,arm_drop,triangle,z,unknown}\r");
            file.Write("@ATTRIBUTE MinX NUMERIC\r");
            file.Write("@ATTRIBUTE MinY NUMERIC\r");
            file.Write("@ATTRIBUTE MinZ NUMERIC\r");
            file.Write("@ATTRIBUTE MaxX NUMERIC\r");
            file.Write("@ATTRIBUTE MaxY NUMERIC\r");
            file.Write("@ATTRIBUTE MaxZ NUMERIC\r");
            file.Write("@ATTRIBUTE MeanX NUMERIC\r");
            file.Write("@ATTRIBUTE MeanY NUMERIC\r");
            file.Write("@ATTRIBUTE MeanZ NUMERIC\r");
            file.Write("@ATTRIBUTE MedianX NUMERIC\r");
            file.Write("@ATTRIBUTE MedianY NUMERIC\r");
            file.Write("@ATTRIBUTE MedianZ NUMERIC\r");
            file.Write("@ATTRIBUTE DiagonalLength NUMERIC\r");
            file.Write("@ATTRIBUTE StartAngleSinXY NUMERIC\r");
            file.Write("@ATTRIBUTE StartAngleCosXY NUMERIC\r");
            file.Write("@ATTRIBUTE StartAngleSinXZ NUMERIC\r");
            file.Write("@ATTRIBUTE FirstLastAngleSinXY NUMERIC\r");
            file.Write("@ATTRIBUTE FirstLastAngleCosXY NUMERIC\r");
            file.Write("@ATTRIBUTE FirstLastAngleSinXZ NUMERIC\r");
            file.Write("@ATTRIBUTE TotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE TotalAngleXZ NUMERIC\r");
            file.Write("@ATTRIBUTE AbsoluteTotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE AbsoluteTotalAngleXZ NUMERIC\r");
            file.Write("@ATTRIBUTE SquaredTotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE SquaredTotalAngleXZ NUMERIC\r");
            file.Write("@ATTRIBUTE FirstLastDistance NUMERIC\r");
            file.Write("@ATTRIBUTE TotalDistance NUMERIC\r");
            file.Write("@ATTRIBUTE MaxAccelerationSquared NUMERIC\r");
            file.Write("@ATTRIBUTE ZeroCrossingsX NUMERIC\r");
            file.Write("@ATTRIBUTE ZeroCrossingsY NUMERIC\r");
            file.Write("@ATTRIBUTE ZeroCrossingsZ NUMERIC\r");
            /*file.Write("@ATTRIBUTE Velocity_MinX NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MinY NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MinZ NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MaxX NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MaxY NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MaxZ NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MeanX NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MeanY NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MeanZ NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MedianX NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MedianY NUMERIC\r");
            file.Write("@ATTRIBUTE Velocity_MedianZ NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MinX NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MinY NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MinZ NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MaxX NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MaxY NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MaxZ NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MeanX NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MeanY NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MeanZ NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MedianX NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MedianY NUMERIC\r");
            file.Write("@ATTRIBUTE InverseVelocity_MedianZ NUMERIC\r");*/
            file.Write("@ATTRIBUTE Stroke_MinX NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MinY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MinZ NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MaxX NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MaxY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MaxZ NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MeanX NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MeanY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MeanZ NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MedianX NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MedianY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MedianZ NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_DiagonalLength NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_StartAngleSinXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_StartAngleCosXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_StartAngleSinXZ NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_FirstLastAngleSinXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_FirstLastAngleCosXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_FirstLastAngleSinXZ NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_TotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_TotalAngleXZ NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_AbsoluteTotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_AbsoluteTotalAngleXZ NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_SquaredTotalAngleXY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_SquaredTotalAngleXZ NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_FirstLastDistance NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_TotalDistance NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_MaxAccelerationSquared NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_ZeroCrossingsX NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_ZeroCrossingsY NUMERIC\r");
            file.Write("@ATTRIBUTE Stroke_ZeroCrossingsZ NUMERIC\r");            
            file.Write("@ATTRIBUTE Duration NUMERIC\r");


            //write out the feature vector for this sample
            file.Write("\r@DATA\r");
            foreach (GestureSample gs in samples)
            {
                string datum = gs.Gesture.ToString();
                for (int i = 0; i < gs.Features.Length; ++i)
                    datum += "," + gs.Features[i];
                file.Write(datum + "\r");
            }

            file.Flush();
            file.Close();
        }
    }
}
