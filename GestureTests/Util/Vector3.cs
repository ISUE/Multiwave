using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

/*
 
Author: Salman Cheema
University of Central Florida
 
Email: salmanc@cs.ucf.edu
 
Released as part of the 3D Gesture Database analysed in
 
"Salman Cheema, Michael Hoffman, Joseph J. LaViola Jr., 3D Gesture classification with linear acceleration and angular velocity 
sensing devices for video games, Entertainment Computing, Volume 4, Issue 1, February 2013, Pages 11-24, ISSN 1875-9521, 10.1016/j.entcom.2012.09.002"
 
*/

namespace GestureTests.Util
{
    /// <summary>
    /// A 3D vector.
    /// </summary>
    public class Vector3
    {
        public float X { get; set; }
        public float Y { get; set; }
        public float Z { get; set; }

        public Vector3()
        {
            X = Y = Z = 0.0f;
        }

        public Vector3(float x, float y, float z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public static Vector3 operator +(Vector3 v1, Vector3 v2)
        {
            return new Vector3(v1.X + v2.X, v1.Y + v2.Y, v1.Z + v2.Z);
        }

        public static Vector3 operator -(Vector3 v1, Vector3 v2)
        {
            return new Vector3(v1.X - v2.X, v1.Y - v2.Y, v1.Z - v2.Z);
        }

        public float LengthSquared()
        {
            return X * X + Y * Y + Z * Z;
        }


        public static Vector3 operator /(Vector3 v1, float f)
        {
            return new Vector3(v1.X / f, v1.Y / f, v1.Z / f);
        }
        public static Vector3 operator *(Vector3 v1, float f)
        {
            return new Vector3(v1.X * f, v1.Y * f, v1.Z * f);
        }

        public static float operator *(Vector3 v1, Vector3 v2)
        {
            return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
        }

        public float Magnitude()
        {
            return (float)Math.Sqrt(LengthSquared());
        }
        public Vector3 Normalize()
        {
            return this / Magnitude();
        }
    }
}
