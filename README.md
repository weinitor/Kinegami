# Kinegami
Given a D-H representation of a kinematic chain robot, the program generates a crease pattern that folds into a kinematically equivalent robot.

Schematic of the numbering of links and joints in a kinematic chain mechanism:
![DH](https://user-images.githubusercontent.com/50150425/129768985-ef83cde7-f979-4d99-b40f-8e7e3fab6ecf.png?v=4&s=10)
<!-- <img src="https://user-images.githubusercontent.com/50150425/129768985-ef83cde7-f979-4d99-b40f-8e7e3fab6ecf.png" width="400" height="200"> -->

The folded state of the origami module, its spacial operator representation, and its crease pattern: (A) the origami prism tube, (B) the flange, (C) the elbow fitting, (D) the twist fitting, (E) the revolute joint, and (E) the prismatic joint.
![OrigamiModule](https://user-images.githubusercontent.com/50150425/129769452-96a78ed9-5819-4436-8d91-4290fb1ca38f.png)


# User Guide:
Run scripts SphericalKinegami.m, UniversalKinegami.m, PlanarKinegami.m, or CylindricalKinegami.m and change parameters if desired. For more comprehensive understanding, reference supporting functions.

# Updates:
7/5/2021:
Edited JointAssignment.m to include correct value of rs for Prismatic Joints.
Edited Kinegami.m to support plotting for Proximal and Distal Frames. Added new function frameplot.m for frame plotting. Changed manner in which figures are closed in papercut files.

# References
Run any files with Kinegami (preferably Universal and Spherical)
Fixing solvedubins3d.m, specifically looking at lines 40-53 (theta calculations)
