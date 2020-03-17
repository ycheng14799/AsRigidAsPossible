import java.awt.*;

public class MyGraphics {
	public Color paintColor, triColor, bkColor, constraintColor, fitColor; 
	public double scale; 
	public int pointRadius; 
	public MyPolygon thePoly; 

	public MyGraphics(Color paint, Color tri, Color bk, Color constraint, Color fit) {
		paintColor = paint; 
		triColor = tri; 
		bkColor = bk;
		constraintColor = constraint; 
		fitColor = fit;
		scale = 1.0D; 
		pointRadius = 4; 
		thePoly = new MyPolygon(); 
	}

	public void clear (Component component, Color paramColor) {
		thePoly = new MyPolygon();
		component.repaint();
	}

	public void redrawPolygon(Graphics g) {
		for(int i=0;i<thePoly.npoints;i++) {
			drawPoint(g, thePoly.xpoints[i], thePoly.ypoints[i]);
		}
		for(int i=0;i<thePoly.npoints-1;i++) {
			g.drawLine(thePoly.xpoints[i], 
						thePoly.ypoints[i],
						thePoly.xpoints[i+1], 
						thePoly.ypoints[i+1]);
		}
		g.drawLine(thePoly.xpoints[thePoly.npoints-1], 
						thePoly.ypoints[thePoly.npoints-1],
						thePoly.xpoints[0], 
						thePoly.ypoints[0]);
	}

	public void triangulate(Graphics g) {
		thePoly.triangulate();
		g.setColor(triColor);
		for(int i=0;i<thePoly.numTriangles;i++) {
			for (int j=0; j <thePoly.triangles[i].npoints; j++) {
				drawPoint(g, thePoly.triangles[i].xpoints[j], thePoly.triangles[i].ypoints[j]);
			}
			for (int j=0; j <thePoly.triangles[i].npoints-1; j++) {
				g.drawLine(thePoly.triangles[i].xpoints[j], 
						thePoly.triangles[i].ypoints[j],
						thePoly.triangles[i].xpoints[j+1], 
						thePoly.triangles[i].ypoints[j+1]);
			}
			g.drawLine(thePoly.triangles[i].xpoints[2], 
				thePoly.triangles[i].ypoints[2],
				thePoly.triangles[i].xpoints[0], thePoly.triangles[i].ypoints[0]);
		}
		g.setColor(paintColor); 
		redrawPolygon(g);
	}

	public void initializeMesh() {
		thePoly.initializeMesh();
	}

	public void addVertex(Graphics g, int x, int y) {
		if (thePoly.npoints > 2) {
			g.setColor(bkColor);
			g.drawLine(thePoly.xpoints[thePoly.npoints-1],
				thePoly.ypoints[thePoly.npoints-1],
				thePoly.xpoints[0],
				thePoly.ypoints[0]);
		}
		g.setColor(paintColor);
		if (thePoly.npoints>=2) {
			g.drawLine(x,y,
				thePoly.xpoints[0],
				thePoly.ypoints[0]);
		} 
		if (thePoly.npoints>0) {
			g.drawLine(thePoly.xpoints[thePoly.npoints-1],
				   thePoly.ypoints[thePoly.npoints-1],
				   x,y);
		}
		thePoly.addPoint(x,y);
		drawPoint(g, x, y);
	}
	
	// Set vertex constraint
	public void addConstraint(Graphics g, int x, int y) {
		int[] constrainedPoint = thePoly.addConstraint(x, y); 
		
		if(constrainedPoint[0] == 1) {
			drawPoint(g, constrainedPoint[1], constrainedPoint[2]);
		}
	}

	
	// Compute components for shape manipulation 
	public void preComputeForManipulation() {
		thePoly.calcStepOneMatrix(); 
		thePoly.calcFInvC(); 
	}

	
	// Get active point for shape manipulation
	public int[] setConstraintActive(int x, int y) {
		return thePoly.setConstraintActive(x, y);
	}

	// Draw updated polygon 
	public void drawUpdatedPoly(Graphics g) {
		g.setColor(paintColor);
		for(int i=0;i<thePoly.updatedPoly.npoints;i++) {
			drawPoint(g, thePoly.updatedPoly.xpoints[i], thePoly.updatedPoly.ypoints[i]);
		}
		for(int i=0;i<thePoly.updatedPoly.npoints-1;i++) {
			g.drawLine(thePoly.updatedPoly.xpoints[i], 
						thePoly.updatedPoly.ypoints[i],
						thePoly.updatedPoly.xpoints[i+1], 
						thePoly.updatedPoly.ypoints[i+1]);
		}
		g.drawLine(thePoly.updatedPoly.xpoints[thePoly.updatedPoly.npoints-1], 
						thePoly.updatedPoly.ypoints[thePoly.updatedPoly.npoints-1],
						thePoly.updatedPoly.xpoints[0], 
						thePoly.updatedPoly.ypoints[0]);
	}
	
	// Draw updated triangulation 
	public void drawUpdatedTriangulation(Graphics g) {
		g.setColor(triColor);
		for(int i=0;i<thePoly.numTriangles;i++) {
			for (int j=0; j <thePoly.updatedTri[i].npoints; j++) {
				drawPoint(g, thePoly.updatedTri[i].xpoints[j], thePoly.updatedTri[i].ypoints[j]);
			}
			for (int j=0; j <thePoly.updatedTri[i].npoints-1; j++) {
				g.drawLine(thePoly.updatedTri[i].xpoints[j], 
						thePoly.updatedTri[i].ypoints[j],
						thePoly.updatedTri[i].xpoints[j+1], 
						thePoly.updatedTri[i].ypoints[j+1]);
			}
			g.drawLine(thePoly.updatedTri[i].xpoints[2], 
				thePoly.updatedTri[i].ypoints[2],
				thePoly.updatedTri[i].xpoints[0], thePoly.updatedTri[i].ypoints[0]);
		}
	}
	
	// Draw final 
	public void drawFinalPoly(Graphics g) {
		g.setColor(paintColor);
		for(int i=0;i<thePoly.finalPoly.npoints;i++) {
			drawPoint(g, thePoly.finalPoly.xpoints[i], thePoly.finalPoly.ypoints[i]);
		}
		for(int i=0;i<thePoly.finalPoly.npoints-1;i++) {
			g.drawLine(thePoly.finalPoly.xpoints[i], 
						thePoly.finalPoly.ypoints[i],
						thePoly.finalPoly.xpoints[i+1], 
						thePoly.finalPoly.ypoints[i+1]);
		}
		g.drawLine(thePoly.finalPoly.xpoints[thePoly.finalPoly.npoints-1], 
						thePoly.finalPoly.ypoints[thePoly.finalPoly.npoints-1],
						thePoly.finalPoly.xpoints[0], 
						thePoly.finalPoly.ypoints[0]);
	}
	
	// Draw updated triangulation 
	public void drawFinalTriangulation(Graphics g) {
		g.setColor(triColor);
		for(int i=0;i<thePoly.numTriangles;i++) {
			for (int j=0; j <thePoly.finalTri[i].npoints; j++) {
				drawPoint(g, thePoly.finalTri[i].xpoints[j], thePoly.finalTri[i].ypoints[j]);
			}
			for (int j=0; j <thePoly.finalTri[i].npoints-1; j++) {
				g.drawLine(thePoly.finalTri[i].xpoints[j], 
						thePoly.finalTri[i].ypoints[j],
						thePoly.finalTri[i].xpoints[j+1], 
						thePoly.finalTri[i].ypoints[j+1]);
			}
			g.drawLine(thePoly.finalTri[i].xpoints[2], 
				thePoly.finalTri[i].ypoints[2],
				thePoly.finalTri[i].xpoints[0], thePoly.finalTri[i].ypoints[0]);
		}
	}

	// redraw Constraints 
	public void redrawConstaints(Graphics g) {
		g.setColor(constraintColor);
		int numConstraints = thePoly.constrainedIdx.size(); 
		for(int i=0; i<numConstraints; i++) {
			drawPoint(g, thePoly.deformedVertices[thePoly.constrainedIdx.get(i)][0],
				thePoly.deformedVertices[thePoly.constrainedIdx.get(i)][1]);
		}
	}
	
	// Manipulate shape 
	public void stepOne(int idx, int newX, int newY) {
		thePoly.stepOne(idx, newX, newY);
	}
	public void drawStepOne(Graphics g) {
		if(thePoly.isUpdated) {
			drawUpdatedTriangulation(g);
			drawUpdatedPoly(g);
			redrawConstaints(g);
		}
	}

	// Fit triangle 
	public void stepTwoOne() {
		thePoly.stepTwoOne();
	}
	public void drawStepTwoOne(Graphics g) {
		if(thePoly.isFit) {
			g.setColor(fitColor);
			for(int i=0; i<thePoly.numTriangles; i++) {
				for (int j=0; j <thePoly.fitTri[i].npoints-1; j++) {
					g.drawLine(thePoly.fitTri[i].xpoints[j], 
							thePoly.fitTri[i].ypoints[j],
							thePoly.fitTri[i].xpoints[j+1], 
							thePoly.fitTri[i].ypoints[j+1]);
				}
				g.drawLine(thePoly.fitTri[i].xpoints[2], 
					thePoly.fitTri[i].ypoints[2],
					thePoly.fitTri[i].xpoints[0], thePoly.fitTri[i].ypoints[0]);
			}
		}
	}
	
	// Final step 
	public void stepTwoTwoSimple() {
		thePoly.stepTwoTwoSimple();
	}
	public void drawStepTwoTwoSimple(Graphics g) {
		if(thePoly.isStepTwoTwo) {
			g.setColor(triColor);
			drawFinalTriangulation(g);
			drawFinalPoly(g);
			redrawConstaints(g);
		}
	}

	// Drawing a point 
	public void drawPoint(Graphics g, int x, int y) {
		g.fillOval(x - pointRadius, 
			y - pointRadius, 
			pointRadius + pointRadius,
			pointRadius + pointRadius); 
	}
}