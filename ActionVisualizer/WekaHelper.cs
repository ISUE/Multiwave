using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Media;
using System.Windows.Input;
using System.Windows;
using System.IO;
using GestureTests.Data;
using GestureTests.Gesture;
using GestureTests.Util;
using GestureTests.Types;
using GestureTests.Experiment;
using GestureTests;
using weka.core;
using weka.classifiers;
using System.Windows.Media.Media3D;


namespace ActionVisualizer
{
    class WekaHelper
    {
        static GestureSample GS;
        static Classifier cls;
        static Classifier cls3D;
        //static ExperimentControl EC;

        public static void initialize()
        {
            //read classifier            
            cls = (Classifier)weka.core.SerializationHelper.read(Config.Classifier2D);
            cls3D = (Classifier)weka.core.SerializationHelper.read(Config.Classifier6D);
            //EC = new ExperimentControl();
            //EC.Initialize();
        }

        public static string Classify(bool useRubine, float duration, bool righthandedness, List<float> SpeakerAngles, PointCollection pointHist, StylusPointCollection S, List<List<int>> hist, List<List<int>> ihist)
        {
            List<Vector2> InterpretedPoints = new List<Vector2>();
            List<Vector2> StylusPoints = new List<Vector2>();
            List<Vector2> VelocityHistory = new List<Vector2>();
            List<Vector2> InverseVelocityHistory = new List<Vector2>();
            foreach(Point P in pointHist)
                InterpretedPoints.Add(new Vector2((float)P.X,(float)P.Y));
            foreach(StylusPoint P in S)
                StylusPoints.Add(new Vector2((float)P.X,(float)P.Y));
            for (int i = 0; i < hist[0].Count; i++)
            {
                VelocityHistory.Add(new Vector2(hist[0][i], hist[1][i]));
                InverseVelocityHistory.Add(new Vector2(ihist[0][i], ihist[1][i]));
            }
            GS = new GestureSample(GestureTests.Types.GestureType.unknown, righthandedness,duration,SpeakerAngles,InterpretedPoints,StylusPoints,VelocityHistory,InverseVelocityHistory);
            GS.ComputeFeatures(GestureFeatures.PointsStroke);

            //if (useRubine)
            //    return EC.Recognizer.Classify(GS).ToString();
            WriteARFF();

            Instances test = new Instances(new java.io.FileReader("outfile.arff"));
            //Instances test = new Instances(new java.io.FileReader("../../../data_arff/u001-training.arff"));
            test.setClassIndex(0);

            double clsLabel = cls.classifyInstance(test.instance(0));
            test.instance(0).setClassValue(clsLabel);
            //Console.WriteLine(clsLabel);

            return ((GestureType2D)((int)clsLabel+1)).ToString();
        }

        public static string Classify3D(bool useRubine, float duration, bool righthandedness, List<float> SpeakerAngles, List<float> SpeakerElevations, Point3DCollection point3DHist, List<List<int>> hist, List<List<int>> ihist)
        {
            List<Vector3> InterpretedPoints = new List<Vector3>();
            List<Vector3> StylusPoints = new List<Vector3>();
            List<Vector3> VelocityHistory = new List<Vector3>();
            List<Vector3> InverseVelocityHistory = new List<Vector3>();
            foreach (Point3D P in point3DHist)
                InterpretedPoints.Add(new Vector3((float)P.X, (float)P.Y, (float)P.Z));

            Vector3 origin = new Vector3(0, 0, 0);
            StylusPoints.Add(origin);
            for (int i = 0; i < InterpretedPoints.Count(); i++)
            {
                origin.X += InterpretedPoints[i].X;
                origin.Y += InterpretedPoints[i].Y;
                origin.Z += InterpretedPoints[i].Z;
                StylusPoints.Add(origin);
            }

            for (int i = 0; i < hist[0].Count; i++)
            {
                VelocityHistory.Add(new Vector3(hist[0][i], hist[1][i], hist[2][i]));
                InverseVelocityHistory.Add(new Vector3(ihist[0][i], ihist[1][i], ihist[2][i]));
            }
            //GestureSample(GestureType gesture, bool rightHanded, float duration, List<float> angles, List<float> elevations, List<Vector3> interpretedPoints, List<Vector3> strokePoints, List<Vector3> velocities, List<Vector3> inverseVelocities)
            GS = new GestureSample(GestureTests.Types.GestureType.unknown, righthandedness, duration, SpeakerAngles, SpeakerElevations, InterpretedPoints, StylusPoints, VelocityHistory, InverseVelocityHistory);
            GS.ComputeFeatures(GestureFeatures.PointsStroke);

            //if (useRubine)
            //    return EC.Recognizer.Classify(GS).ToString();
            WriteARFF3D();

            Instances test = new Instances(new java.io.FileReader("outfile.arff"));
            //Instances test = new Instances(new java.io.FileReader("../../../data_arff/u001-training.arff"));
            test.setClassIndex(0);

            double clsLabel = cls3D.classifyInstance(test.instance(0));
            test.instance(0).setClassValue(clsLabel);
            //Console.WriteLine(clsLabel);
            //Console.WriteLine(((GestureType)((int)clsLabel+1)).ToString());
            
            return ((GestureType)((int)clsLabel + 1)).ToString();
        }

        private static void WriteARFF()
        {
            //create file
            StreamWriter file = File.CreateText("outfile.arff");

            //write arff header
            file.Write("@RELATION " + "TestSample" + "\r\r");

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
            file.Write("@ATTRIBUTE Inverse_MaxAccelerationSquared NUMERIC\r");/*
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

            string datum = GS.Gesture.ToString();
            for (int i = 0; i < GS.Features.Length; ++i)
                datum += "," + GS.Features[i];
            file.Write(datum + "\r");
            file.Flush();
            file.Close();
        }

        private static void WriteARFF3D()
        {
            //create file
            StreamWriter file = File.CreateText("outfile.arff");

            //write arff header
            file.Write("@RELATION " + "TestSample" + "\r\r");

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

            string datum = GS.Gesture.ToString();
            for (int i = 0; i < GS.Features.Length; ++i)
                datum += "," + GS.Features[i];
            file.Write(datum + "\r");

            file.Flush();
            file.Close();
        }
    }
}
