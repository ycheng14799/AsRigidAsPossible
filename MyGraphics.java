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

		// TODO: Case of outer boundary intersecting triangulation
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