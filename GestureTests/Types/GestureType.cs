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


namespace GestureTests.Types
{
    /// <summary>
    /// An enumeration of the different types of supported gestures.
    /// </summary>

    public enum GestureType
    {
        triangle = 1,
        x = 2,
        rectangle = 3,
        circle = 4,
        c = 5,
        check = 6,
        caret = 7,
        zigzag = 8,
        arrow = 9,
        star = 10,
        z = 11,
        w = 12,
        double_arch = 13,
        mu = 14,
        y = 15,
        s = 16,
        p = 17,
        q = 18,        
        unknown = 19,
        side_x = 20,
        left_bracket = 21,
        right_bracket = 22,
    };
    public enum GestureType2D
    {
        triangle = 1,
        x = 2,
        side_x = 3,
        rectangle = 4,
        circle = 5,
        check = 6,
        caret = 7,
        zigzag = 8,
        arrow = 9,
        left_bracket = 10,
        right_bracket = 11,
        star = 12,
        unknown = 13
    };
}
