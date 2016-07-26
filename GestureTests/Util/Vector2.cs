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
    /// A 2D vector.
    /// </summary>
    public class Vector2
    {
        public float X { get; set; }
        public float Y { get; set; }        

        public Vector2()
        {
            X = Y = 0.0f;
        }

        public Vector2(float x, float y)
        {
            X = x;
            Y = y;
        }

        public static Vector2 operator +(Vector2 v1, Vector2 v2)
        {
            return new Vector2(v1.X + v2.X, v1.Y + v2.Y);
        }

        public static Vector2 operator -(Vector2 v1, Vector2 v2)
        {
            return new Vector2(v1.X - v2.X, v1.Y - v2.Y);
        }

        public static Vector2 operator /(Vector2 v1, float f)
        {
            return new Vector2(v1.X/f, v1.Y/f);
        }
        public static Vector2 operator *(Vector2 v1, float f)
        {
            return new Vector2(v1.X * f, v1.Y * f);
        }

        public static float operator *(Vector2 v1, Vector2 v2)
        {
            return v1.X * v2.X + v1.Y * v2.Y;
        }

        public float LengthSquared()
        {
            return X * X + Y * Y;
        }

        public float Magnitude()
        {
            return (float)Math.Sqrt(LengthSquared());
        }
        public Vector2 Normalize()
        {
            return this/Magnitude();
        }
    }
}
