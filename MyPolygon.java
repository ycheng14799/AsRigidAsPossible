import java.awt.*;
import Jama.*; 

public class MyPolygon extends Polygon {
	public Polygon[] triangles;
	public int[] tempX = new int[3];
	public int[] tempY = new int[3];
	public int numTriangles; 
	public boolean triangulated = false; 
	public int[][] gComponentsSingle = new int[6][6];
	//public int[] constraintX = new int[3];
	//public int[] constraintY = new int[3]; 

	private void resizeTriangles(int num) {
		assert(tempX.length == tempY.length); 
		while(3*num > tempX.length) {
			// Double array capacity 
			int[] copyX = new int[tempX.length*2]; 
			int[] copyY = new int[tempX.length*2];
			// Copy data across 
			System.arraycopy(tempX, 0, copyX, 0, tempX.length);
			System.arraycopy(tempY, 0, copyY, 0, tempX.length);
			// Reset references 
			tempX = copyX; 
			tempY = copyY; 
		}
		assert(tempX.length == tempY.length);
	}

	private void resizeConstraints(int num) {
	/*
		assert(constraintX.length == constraintY.length); 
		while(num > constraintX.length) {
			// Double array capacity
			int[] copyX = new int[constraintX.length*2]; 
			int[] copyY = new int[constraintX.length*2];
			// Copy data across 
			System.arraycopy(constraintX, 0, copyX, 0, tempX.length);
			System.arraycopy(constraintY, 0, copyY, 0, tempX.length);
			// Reset references 
			constraintX = copyX;
			constraintY = copyY;
		}
		assert(constraintX.length == constraintY.length);
	*/
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
		double cx, cy;
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
							// Exclude triangles entirely outside 
							cx = (xpoints[i] + xpoints[j] + xpoints[k]) / 3.0;
							cy = (ypoints[i] + ypoints[j] + ypoints[k]) / 3.0;
							if(contains(cx,cy)) {
								numTriangles++; 
								resizeTriangles(numTriangles);
								tempX[3*(numTriangles-1)] = xpoints[i];
								tempX[3*(numTriangles-1) + 1] = xpoints[j];
								tempX[3*(numTriangles-1) + 2] = xpoints[k];
								tempY[3*(numTriangles-1)] = ypoints[i];
								tempY[3*(numTriangles-1) + 1] = ypoints[j];
								tempY[3*(numTriangles-1) + 2] = ypoints[k];
								tempY[3*(numTriangles-1) + 2] = ypoints[k];
							}
						}
					}
				}
			}
		}
		triangles = new Polygon[numTriangles];
		int[] newPolyX = new int[3]; 
		int[] newPolyY = new int[3];
		for(int i=0; i<numTriangles; i++) {
			for(int j=0; j<3; j++) {
				newPolyX[j] = tempX[3*i + j];
				newPolyY[j] = tempY[3*i + j]; 
			}
			triangles[i] = new Polygon(newPolyX, newPolyY, 3);
		}
		triangulated = true; 
	}

	// Precompute components for scale-free manipulation
	// Step One in Igarashi et. al's paper
	// Helper method for calculating 
	// matrix components for a single triangle
	// @params: v0x, v0y, v1x, v1y, v2x, v2y 
	// Output: integer array of matrix components 
	public void calcComponentsSingle(int v0x, int v0y, int v1x, int v1y, int v2x, int v2y) {

		// find x01, y01
		int x01 = (v2x-v0x)*(v1x-v0x) + (v2y-v0y)*(v1y-v0y);
		int y01 = (v2x-v0x)*(v1y-v0y) + (v2y-v0y)*(v0x-v1x);

		// Populate gComponentsSingle[] array 
		gComponentsSingle[0][0] = x01*x01 - 2*x01 + y01*y01 + 1; 
		gComponentsSingle[2][0] = - x01*x01 + 2*x01 - 2*y01*y01;
		gComponentsSingle[3][0] = 2*y01; 
		gComponentsSingle[4][0] = 2*x01 - 2; 
		gComponentsSingle[5][0] = - 2*y01; 
		gComponentsSingle[1][1] = x01*x01 - 2*x01 + y01*y01 + 1; 
		gComponentsSingle[2][1] = - 2*y01; 
		gComponentsSingle[3][1] = - 2*y01*y01 + 2*x01 - 2*x01*x01; 
		gComponentsSingle[4][1] = 2*y01; 
		gComponentsSingle[5][1] = 2*x01 - 2; 
		gComponentsSingle[2][2] = x01*x01 + y01*y01; 
		gComponentsSingle[4][2] = - 2*x01; 
		gComponentsSingle[5][2] = 2*y01; 
		gComponentsSingle[3][3] = x01*x01 + y01*y01; 
		gComponentsSingle[4][3] = - 2*y01; 
		gComponentsSingle[5][3] = - 2*x01; 
		gComponentsSingle[4][4] = 1; 
		gComponentsSingle[5][5] = 1;
	}

	public void calcGMatrix() { 
		
	}


}