+++
date = "2017-05-05T22:55:50+05:30"
title = "Meshing Techniques - Cambered Airfoil"
draft = true

+++

Meshing is a crucial process in obtaining accurate results for various simulations across different fields. In computational fluid dynamics, various meshing techniques are used in grid generation for 2D analyses of airfoils. Some nice runthroughs exist on YouTube, but they mostly deal with symmetric airfoils such as the *beloved* NACA 0012. I haven't found any nice documentation that covers standard techniques for different kinds of airfoils, so I thought I'd start with some.

---

**Airfoil Selection:** A new airfoil I've been researching is the Selig S7075 for its good gliding/soaring characteristics and high efficiency at low Reynolds numbers. It's a cambered airfoil with sharp trailing edge, and a maximum thickness-to-chord ratio of 9%. XFOIL analyses predict accurate results in agreement with Dr. Selig's [A Summary of Low-Speed Airfoil Data, Vol. 2](https://m-selig.ae.illinois.edu/uiuc_lsat/Low-Speed-Airfoil-Data-V2.pdf), whose experiments were conducted using UIUC's wind tunnel. These two sources of data should be adequate for verification and validation. Let's see what CFD results have to offer.

---

**Meshing Software:** ANSYS Workbench Mesh is terrible at meshing cambered airfoils, and flow analyses from these meshes do not provide accurate results.
ICEM CFD is a better option because it is far more customisable and generates good meshes quickly.

The two types of meshes I'll be analysing are an O-grid and a C-grid. An O-grid setup is when the body of analysis is enclosed within a curvilinear grid in which completely closed contours are created at offsets from the body. A C-grid setup, by example for an airfoil, is when the surrounding grid forms a 'C' till the trailing edge.  Forums recommend that a C-grid mesh provides better results for airfoils with sharp trailing edges. This is because the grid alignment is along the airfoil shape (hence the flow streamlines) to capture the flow results more accurately.

The blocking process in ICEM CFD used the same ratios for element sizing in both meshing techniques.

---
**C-Grid:**

<div style="text-align:center"><img src="../Meshing Techniques - Cambered Airfoil/C-Grid.png" width="95%"></div>
<div style="text-align:center"><img src="../Meshing Techniques - Cambered Airfoil/C-GridBound.png" width="95%"></div>
<div style="text-align:center"><img src="../Meshing Techniques - Cambered Airfoil/C-GridLE.png" width="95%"></div>

---
**O-Grid:**

<div style="text-align:center"><img src="../Meshing Techniques - Cambered Airfoil/O-Grid.png" width="95%"></div>
<div style="text-align:center"><img src="../Meshing Techniques - Cambered Airfoil/O-GridBound.png" width="95%"></div>
<div style="text-align:center"><img src="../Meshing Techniques - Cambered Airfoil/O-GridLE.png" width="95%"></div>
