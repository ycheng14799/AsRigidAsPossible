import java.awt.*;

public class MyGraphics {
	public Color paintColor, triColor, bkColor; 
	public double scale; 
	public int pointRadius; 
	public MyPolygon thePoly; 

	public MyGraphics(Color paint, Color tri, Color bk) {
		paintColor = paint; 
		triColor = tri; 
		bkColor = bk;
		scale = 1.0D; 
		pointRadius = 4; 
		thePoly = new MyPolygon(); 
	}

	public void clear (Component component, Color paramColor) {
		thePoly = new MyPolygon();
		component.repaint();
	}

	public void triangulate(Graphics g) {
		thePoly.triangulate();
		/*
		for(int i=0;i<thePoly.numTriangles;i++) {
			g.setColor(triColor);
			for (int j=0; j <thePoly.triangles[i].npoints; j++) {
				drawPoint(g, thePoly.triangles[i].xpoints[j], thePoly.triangles[i].ypoints[j]);
			}
			//g.drawLine(thePoly.triangles[i].xpoints[thePoly.numTriangles-1], 
			//	thePoly.triangles[i].ypoints[thePoly.numTriangles-1],
			//	thePoly.triangles[i].xpoints[0], thePoly.triangles[i].ypoints[0]);
		}
		*/
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

	public void drawPoint(Graphics g, int x, int y) {
		g.fillOval(x - pointRadius, 
			y - pointRadius, 
			pointRadius + pointRadius,
			pointRadius + pointRadius); 
	}
}