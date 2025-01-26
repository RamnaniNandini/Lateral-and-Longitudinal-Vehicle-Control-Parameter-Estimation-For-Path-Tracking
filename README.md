# Lateral-and-Longitudinal-Vehicle-Control-Parameter-Estimation-For-Path-Tracking
Path tracking is a crucial aspect of autonomous vehicle control, enabling vehicles to follow
predetermined trajectories accurately and safely. This research aims to improve how
these cars control their sideways movement (lateral) and forward speed (longitudinal).
By combining advanced control techniques with dynamic vehicle modeling, we aim to
enhance the performance and reliability of autonomous navigation systems.

The primary objective of this project was to optimize various parameters—such as heading
error, lateral error, steering angle, and yaw rate—to achieve accurate velocity control and
precise path following. The final model employed a Linear Quadratic Regulator (LQR)
approach to optimize these parameters and fit the path. This method utilized a grid
optimization technique to determine appropriate values for Q and R matrices, optimizing
the controller gain (K) for both case of constant and varying velocity. Before settling on
this approach, two other methods were thoroughly analyzed: tracking error computation
using discrete waypoint representation and minimizing the distance between actual and
model coordinates.

The project also focused on determining maximum deceleration for various curves and
analyzing them using parameters such as radius of curvature, curve length, maximum
velocity within the curve, and curve direction. To gain deeper insights into vehicle
dynamics, data analysis techniques such as correlation matrices and multiple linear
regression were used.


