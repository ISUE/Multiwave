using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GestureTests.Gesture;
using GestureTests.Types;

/*
 
Author: Salman Cheema
University of Central Florida
 
Email: salmanc@cs.ucf.edu
 
Released as part of the 3D Gesture Database analysed in
 
"Salman Cheema, Michael Hoffman, Joseph J. LaViola Jr., 3D Gesture classification with linear acceleration and angular velocity 
sensing devices for video games, Entertainment Computing, Volume 4, Issue 1, February 2013, Pages 11-24, ISSN 1875-9521, 10.1016/j.entcom.2012.09.002"
 
*/

namespace GestureTests
{
    /// <summary>
    /// configuration settings for training/recognition/experiment setup.
    /// </summary>
    public class Config
    {
        /// <summary>
        /// Source location for dataset to be used for experiment.
        /// </summary>
        public static string DataPath = @"..\..\..\data\";

        /// <summary>
        /// Target location where experiment results will be stored as text files.
        /// </summary>
        public static string ResultsPath = @"..\..\..\..\results\";

        /// <summary>
        /// Target location where the dataset will be exported as '.arff' files to be used for weka.
        /// </summary>
        public static string WekaOutputPath = @"..\..\..\data_arff\";

        /// <summary>
        /// Number of training samples per gesture to be used for the experiment. 
        /// In 'UserDependent' mode, this can range between 0-25, as the number of available samples is 25 samples per user per gesture.
        /// In 'UserIndependent' mode, this can range between 0-625, as the number of available samples is 625 samples for each gesture (from all 25 users).
        /// </summary>
        public static int NumTrainingSamples = 0;

        /// <summary>
        /// Number of times to run the experiment with random sub-selections of the training data.
        /// </summary>
        public static int NumExperiments = 20;

        public static bool ReportUserSpecificAccuracyResults = false;
        public static bool ReportIndividualUsersResults = false;
        
        public static bool Use3DMode = false;

        public static String Classifier2D = "../../../data_arff/RandomForestu010.model";
        public static String Classifier6D = "../../../data_arff/RandomForest6Du010.model";

        public static GestureFeatures FeaturesToUse = GestureFeatures.PointsStroke;

        /// <summary>
        /// Can be used to select a subset of available gestures for the experiment.
        /// </summary>
        public static List<GestureType> GesturesToUse = new List<GestureType>();

        /// <summary>
        /// The number of features are dependent on whether only Wiimote features are to be used (29 features) or whether both Wiimote and MotionPlus (29 + 12) features are to be used.
        /// </summary>
        public static int NumFeatures
        {
            get
            {
                if (Use3DMode)
                {
                    switch (FeaturesToUse)
                    {
                        case GestureFeatures.Points:
                            return 29+3;
                        case GestureFeatures.PointsStroke:
                            return 57+6;
                        case GestureFeatures.PointsStrokeInverse:
                            return 69+6;
                        case GestureFeatures.PointsVelocities:
                            return 41+3;
                        case GestureFeatures.PointsVelocitiesInverseVelocities:
                            return 53+3;
                    }
                    return 81+6;
                }
                else
                {
                    switch (FeaturesToUse)
                    {
                        case GestureFeatures.Points:
                            return 20+2;
                        case GestureFeatures.PointsStroke:
                            return 39+2;
                        case GestureFeatures.PointsStrokeInverse:
                            return 47+2;
                        case GestureFeatures.PointsVelocities:
                            return 28+2;
                        case GestureFeatures.PointsVelocitiesInverseVelocities:
                            return 36+2;
                    }
                    return 63+2;
                }
            }
        }

        static Config()
        {
            if (Use3DMode == false)
            {
                GesturesToUse.Add(GestureType.triangle);
                GesturesToUse.Add(GestureType.x);
                GesturesToUse.Add(GestureType.side_x);
                GesturesToUse.Add(GestureType.rectangle);
                GesturesToUse.Add(GestureType.circle);
                GesturesToUse.Add(GestureType.check);
                GesturesToUse.Add(GestureType.caret);
                GesturesToUse.Add(GestureType.zigzag);
                GesturesToUse.Add(GestureType.arrow);
                GesturesToUse.Add(GestureType.x);
                GesturesToUse.Add(GestureType.left_bracket);
                GesturesToUse.Add(GestureType.right_bracket);   
                GesturesToUse.Add(GestureType.star);
                GesturesToUse.Add(GestureType.unknown);
            }
            else
            {
                GesturesToUse.Add(GestureType.triangle);
                GesturesToUse.Add(GestureType.x);
                GesturesToUse.Add(GestureType.side_x);
                GesturesToUse.Add(GestureType.rectangle);
                GesturesToUse.Add(GestureType.circle);
                GesturesToUse.Add(GestureType.check);
                GesturesToUse.Add(GestureType.caret);
                GesturesToUse.Add(GestureType.zigzag);
                GesturesToUse.Add(GestureType.arrow);
                GesturesToUse.Add(GestureType.x);
                GesturesToUse.Add(GestureType.left_bracket);
                GesturesToUse.Add(GestureType.right_bracket);
                GesturesToUse.Add(GestureType.star);
                GesturesToUse.Add(GestureType.unknown);
            }
     
        }
    }
}
