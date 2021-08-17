# Kinegami
Given a D-H representation of a kinematic chain robot, the program generates a crease pattern that folds into a kinematically equivalent robot.

The numbering system of the D-H representations are show as:
![DH](https://user-images.githubusercontent.com/50150425/129768985-ef83cde7-f979-4d99-b40f-8e7e3fab6ecf.png)

The origami modules used in this program are:
![FittingCP2](https://user-images.githubusercontent.com/50150425/129769294-483909e7-8157-471a-84a9-f250a32967cb.png)

# User Guide:
Run scripts SphericalKinegami.m, UniversalKinegami.m, PlanarKinegami.m, or CylindricalKinegami.m and change parameters if desired. For more comprehensive understanding, reference supporting functions.

# Updates:
7/5/2021:
Edited JointAssignment.m to include correct value of rs for Prismatic Joints.
Edited Kinegami.m to support plotting for Proximal and Distal Frames. Added new function frameplot.m for frame plotting. Changed manner in which figures are closed in papercut files.

# References
Run any files with Kinegami (preferably Universal and Spherical)
Fixing solvedubins3d.m, specifically looking at lines 40-53 (theta calculations)
