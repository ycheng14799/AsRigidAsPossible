import java.awt.*;

public class MyPolygon extends Polygon {
	public Polygon[] triangles;
	public int[] tempX = new int[3];
	public int[] tempY = new int[3];
	public int numTriangles; 
	public boolean triangulated; 

	private void resizeTriangles(int num) {
		assert(tempX.length == tempY.length); 
		while(3*num > tempTri.length) {
			// Double array capacity 
			int[] copyX = new int[tempX.length*2]; 
			int[] copyY = new int[tempX.length*2];
			// Copy data across 
			System.arraycopy(tempX, 0, copyX, tempX.length);
			System.arraycopy(tempY, 0, copyY, tempY.length);
			// Reset references 
			tempX = copyX; 
			tempY = copyY; 
		}
		assert(tempX.length == tempY.length);
	}

	public MyPolygon() {}
	public MyPolygon(int x[], int y[], int n) {}

	// O(n^4) implementation 
	public void triangulate() {
		numTriangles = 0; 

		// Compute z = x^2 + y^2
		int[] z = new int[npoints]; 
		for (int i=0; i<npoints; i++) {
			z[i] = (xpoints[i]*xpoints[i]) +
				(ypoints[i]*ypoints[i]);
		}
		// For each triple (i,j,k)
		int xn, yn, zn;
		boolean flag; 
		for (int i=0; i<npoints-2; i++) {
			for (int j=i+1; j<npoints; j++) {
				for (int k=i+1; k<npoints; k++) {
					// Compute normal 
					xn = (ypoints[j]-ypoints[i])*(z[k]-z[i]) -
						(ypoints[k]-ypoints[i])*(z[j]-z[i]);
					yn = (xpoints[k]-xpoints[i])*(z[j]-z[i]) -
						(xpoints[j]-xpoints[i])*(z[k]-z[i]);
					zn = (xpoints[j]-xpoints[i])*(ypoints[k]-ypoints[i]) -
						(xpoints[k]-xpoints[i])*(ypoints[j]-ypoints[i]);
					// Examine faces on the bottom of the paraboloid
					if(zn < 0) {
						flag = true; 
						for(int m=0; m<npoints; m++) {
							flag = flag &&
								((xpoints[m]-xpoints[i])*xn +
								(ypoints[m]-ypoints[i])*yn +
								(z[m]-z[i])*zn <= 0);
						}
						if (flag) {
							numTriangles++; 
							resizeTriangles(numTriangles);
							tempX[3*(numTriangles-1)] = xpoints[i];
							tempX[3*(numTriangles-1) + 1] = xpoints[j];
							tempX[3*(numTriangles-1) + 2] = xpoints[k];
							/*
							numTriangles++;
							resizeTriangles(numTriangles);
							triangles[numTriangles-1] = new Polygon(new int[]{xpoints[i],xpoints[j],xpoints[k]},
								new int[]{ypoints[i],ypoints[j],ypoints[k]}, 3);
							*/
						}
					}
				}
			}
		}
	}
}