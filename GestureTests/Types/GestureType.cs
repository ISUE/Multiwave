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
        swipe_left = 13,
        swipe_right = 14,
        swipe_up = 15,
        swipe_down = 16,
        swipe_front = 17,
        swipe_back = 18,
        tap_left = 19,
        tap_right = 20,
        tap_up = 21,
        tap_down = 22,
        tap_front = 23,
        tap_back = 24,
        scratchout = 25,
        square = 26,
        c = 27,
        two_handed_fb = 28,
        two_handed_lr = 29,
        horizontal_circle = 30,
        vertical_circle = 31,
        spiral = 32,
        arm_lift = 33,
        z = 34,
        unknown = 35
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
