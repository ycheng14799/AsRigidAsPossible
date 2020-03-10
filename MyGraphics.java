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

	public void addVertex(Graphics g, int x, int y) {
		Color color = g.getColor(); 
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
		g.setColor(color);
		thePoly.addPoint(x,y);
		drawPoint(g, x, y);
	}

	public void drawPoint(Graphics g, int x, int y) {
		Color color = g.getColor(); 
		g.setColor(paintColor);
		g.fillOval(x - pointRadius, 
			y - pointRadius, 
			pointRadius + pointRadius,
			pointRadius + pointRadius); 
		g.setColor(color);
	}
}