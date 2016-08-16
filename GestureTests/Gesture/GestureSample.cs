using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
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

namespace GestureTests.Gesture
{
    /// <summary>
    /// Contains the datapoints corresponding to a single instance of a gesture. 
    /// </summary>
    public class GestureSample
    {

        /// <summary>
        /// Type of Gesture.
        /// </summary>
        public GestureType Gesture { get; protected set; }

        /// <summary>
        /// Gesture duration in ms.
        /// </summary>
        public float Duration { get; protected set; }

        /// <summary>
        /// List of speaker angles
        /// </summary>
        public List<float> SpeakerAngles { get; protected set; }

        /// <summary>
        /// List of points determined from instantaneous velocities
        /// </summary>
        public List<Vector2> InterpretedPoints { get; protected set; }

        /// <summary>
        /// List of points determined from instantaneous velocities
        /// </summary>
        public List<Vector2> StrokePoints { get; protected set; }

        /// <summary>
        /// List of raw velocities based on bandwidth
        /// </summary>
        public List<Vector2> Velocities { get; protected set; }

        /// <summary>
        /// List of raw inverse velocities based on bandwidth
        /// </summary>
        public List<Vector2> InverseVelocities { get; protected set; }

        public List<Vector2> InversePoints{ get; protected set; }

        public List<float[]> RawData { get; protected set; }


        public List<float> SpeakerElevations { get; protected set; }
        public List<Vector3> InterpretedPoints3D { get; protected set; }
        public List<Vector3> StrokePoints3D { get; protected set; }
        public List<Vector3> Velocities3D { get; protected set; }
        public List<Vector3> InverseVelocities3D { get; protected set; }
        public List<Vector3> InversePoints3D { get; protected set; }

        /// <summary>
        /// A boolean flag indicating whether this gesture was performed with the Right or Left hand.
        /// </summary>
        public bool RightHanded { get; protected set; }

        /// <summary>
        /// A list of features computed from the data points.
        /// </summary>
        public float[] Features { get; protected set; }

        //sample = new GestureSample(gesture, rightHanded, duration, angles, interpretedPoints, velocities, inverseVelocities);
        public GestureSample(GestureType gType, bool rightHanded, float duration, List<float> angles, List<Vector2> interpretedPoints, List<Vector2> strokePoints, List<Vector2> velocities, List<Vector2> inverseVelocities, List<float[]> raw_data)
        {
            Gesture = gType;
            Duration = duration;
            RightHanded = rightHanded;

            SpeakerAngles = new List<float>();
            SpeakerAngles.AddRange(angles);

            InterpretedPoints = new List<Vector2>();
            InterpretedPoints.AddRange(interpretedPoints);

            StrokePoints = new List<Vector2>();
            StrokePoints.AddRange(strokePoints);

            Velocities = new List<Vector2>();
            Velocities.AddRange(velocities);

            InverseVelocities = new List<Vector2>();
            InverseVelocities.AddRange(inverseVelocities);

            InversePoints = new List<Vector2>();
            foreach (Vector2 velocity_pair in InverseVelocities)
            {
                Vector2 P = new Vector2();
                P.X = (float)(velocity_pair.X * Math.Sin(SpeakerAngles[0] * Math.PI / 180.0));
                P.Y = (float)(velocity_pair.X * Math.Cos(SpeakerAngles[0] * Math.PI / 180.0));
                P.X += (float)(velocity_pair.Y * Math.Sin(SpeakerAngles[1] * Math.PI / 180.0));
                P.Y += (float)(velocity_pair.Y * Math.Cos(SpeakerAngles[1] * Math.PI / 180.0));
                InversePoints.Add(P);
            }

            RawData = new List<float[]>();
            RawData.AddRange(raw_data);
        }

        public GestureSample(GestureType gesture, bool rightHanded, float duration, List<float> angles, List<float> elevations, List<Vector3> interpretedPoints, List<Vector3> strokePoints, List<Vector3> velocities, List<Vector3> inverseVelocities, List<float[]> raw_data)
        {
            // TODO: Complete member initialization
            Gesture = gesture;
            Duration = duration;
            RightHanded = rightHanded;

            SpeakerAngles = new List<float>();
            SpeakerAngles.AddRange(angles);
            
            SpeakerElevations = new List<float>();
            SpeakerElevations.AddRange(elevations);

            InterpretedPoints3D = new List<Vector3>();
            InterpretedPoints3D.AddRange(interpretedPoints);

            StrokePoints3D = new List<Vector3>();
            StrokePoints3D.AddRange(strokePoints);

            Velocities3D = new List<Vector3>();
            Velocities3D.AddRange(velocities);

            InverseVelocities3D = new List<Vector3>();
            InverseVelocities3D.AddRange(inverseVelocities);

            InversePoints3D = new List<Vector3>();
            foreach (Vector3 velocity_tuple in InverseVelocities3D)
            {
                Vector3 P = new Vector3();
                P.X = (float)(velocity_tuple.X * Math.Sin(SpeakerAngles[0] * Math.PI / 180.0) * Math.Cos(SpeakerElevations[0] * Math.PI / 180.0));
                P.Y = (float)(velocity_tuple.X * Math.Cos(SpeakerAngles[0] * Math.PI / 180.0) * Math.Cos(SpeakerElevations[0] * Math.PI / 180.0));
                P.Z = (float)(velocity_tuple.X * Math.Sin(SpeakerElevations[0] * Math.PI / 180.0));
                P.X = (float)(velocity_tuple.Y * Math.Sin(SpeakerAngles[1] * Math.PI / 180.0) * Math.Cos(SpeakerElevations[1] * Math.PI / 180.0));
                P.Y = (float)(velocity_tuple.Y * Math.Cos(SpeakerAngles[1] * Math.PI / 180.0) * Math.Cos(SpeakerElevations[1] * Math.PI / 180.0));
                P.Z = (float)(velocity_tuple.Y * Math.Sin(SpeakerElevations[1] * Math.PI / 180.0));
                P.X = (float)(velocity_tuple.Z * Math.Sin(SpeakerAngles[2] * Math.PI / 180.0) * Math.Cos(SpeakerElevations[2] * Math.PI / 180.0));
                P.Y = (float)(velocity_tuple.Z * Math.Cos(SpeakerAngles[2] * Math.PI / 180.0) * Math.Cos(SpeakerElevations[2] * Math.PI / 180.0));
                P.Z = (float)(velocity_tuple.Z * Math.Sin(SpeakerElevations[2] * Math.PI / 180.0));
                
                InversePoints3D.Add(P);
            }


            RawData = new List<float[]>();
            RawData.AddRange(raw_data);
        }

        /// <summary>
        /// Will recompute the feature vector.
        /// </summary>
        /// <param name="typeOfFeaturesToCompute">indicates the type of features to be computed.</param>
        public void ComputeFeatures(GestureFeatures typeOfFeaturesToCompute)
        {
            if (Config.Use3DMode)
            {
                Compute3DFeatures(typeOfFeaturesToCompute);
                return;
            }

            List<float> allFeatures = new List<float>();
            XYFeatures points = null, inversePoints = null, velocities = null, inverseVelocities = null;
            StrokeFeatures strokes = null;

            switch (typeOfFeaturesToCompute)
            {
                case GestureFeatures.Points:
                    points = new XYFeatures(InterpretedPoints);
                    break;
                case GestureFeatures.PointsStroke:
                    points = new XYFeatures(InterpretedPoints);
                    strokes = new StrokeFeatures(StrokePoints);
                    break;
                case GestureFeatures.PointsStrokeInverse:
                    points = new XYFeatures(InterpretedPoints);
                    strokes = new StrokeFeatures(StrokePoints);
                    inversePoints = new XYFeatures(InterpretedPoints);
                    break;
                case GestureFeatures.PointsVelocities:
                    points = new XYFeatures(InterpretedPoints);
                    velocities = new XYFeatures(Velocities);
                    break;
                case GestureFeatures.PointsVelocitiesInverseVelocities:
                    points = new XYFeatures(InterpretedPoints);
                    velocities = new XYFeatures(Velocities);
                    inverseVelocities = new XYFeatures(InverseVelocities);
                    break;
                case GestureFeatures.PointsStrokesVelocitiesInverseVelocities:
                    points = new XYFeatures(InterpretedPoints);
                    velocities = new XYFeatures(Velocities);
                    inverseVelocities = new XYFeatures(InverseVelocities);
                    strokes = new StrokeFeatures(StrokePoints);
                    break;
            }

            //allFeatures.AddRange(SpeakerAngles);
            allFeatures.AddRange(points.NormalFeatures);
            if (inversePoints != null)
                allFeatures.AddRange(inversePoints.VelocityFeatures);
            if (velocities != null)
                allFeatures.AddRange(velocities.VelocityFeatures);
            if (inverseVelocities != null)
                allFeatures.AddRange(inverseVelocities.VelocityFeatures);
            if (strokes != null)
                allFeatures.AddRange(strokes.NormalFeatures);

            allFeatures.Add(this.Duration);
            Features = allFeatures.ToArray();
        }

        private void Compute3DFeatures(GestureFeatures typeOfFeaturesToCompute)
        {
            List<float> allFeatures = new List<float>();
            XYZFeatures points = null, inversePoints = null, velocities = null, inverseVelocities = null, strokes = null;

            switch (typeOfFeaturesToCompute)
            {
                case GestureFeatures.Points:
                    points = new XYZFeatures(InterpretedPoints3D);
                    break;
                case GestureFeatures.PointsStroke:
                    points = new XYZFeatures(InterpretedPoints3D);
                    strokes = new XYZFeatures(StrokePoints3D);
                    break;
                case GestureFeatures.PointsStrokeInverse:
                    points = new XYZFeatures(InterpretedPoints3D);
                    strokes = new XYZFeatures(StrokePoints3D);
                    inversePoints = new XYZFeatures(InterpretedPoints3D);
                    break;
                case GestureFeatures.PointsVelocities:
                    points = new XYZFeatures(InterpretedPoints3D);
                    velocities = new XYZFeatures(Velocities3D);
                    break;
                case GestureFeatures.PointsVelocitiesInverseVelocities:
                    points = new XYZFeatures(InterpretedPoints3D);
                    velocities = new XYZFeatures(Velocities3D);
                    inverseVelocities = new XYZFeatures(InverseVelocities3D);
                    break;
                case GestureFeatures.PointsStrokesVelocitiesInverseVelocities:
                    points = new XYZFeatures(InterpretedPoints3D);
                    velocities = new XYZFeatures(Velocities3D);
                    inverseVelocities = new XYZFeatures(InverseVelocities3D);
                    strokes = new XYZFeatures(StrokePoints3D);
                    break;
            }

            //allFeatures.AddRange(SpeakerAngles);
            allFeatures.AddRange(points.NormalFeatures);
            if (inversePoints != null)
                allFeatures.AddRange(inversePoints.VelocityFeatures);
            if (velocities != null)
                allFeatures.AddRange(velocities.VelocityFeatures);
            if (inverseVelocities != null)
                allFeatures.AddRange(inverseVelocities.VelocityFeatures);
            if (strokes != null)
                allFeatures.AddRange(strokes.NormalFeatures);

            allFeatures.Add(this.Duration);
            Features = allFeatures.ToArray();
        }
    }
}
