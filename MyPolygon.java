import java.awt.*;
import Jama.*; 
import java.util.ArrayList;

public class MyPolygon extends Polygon {
	public Polygon[] triangles;
	public int[] tempX = new int[3];
	public int[] tempY = new int[3];
	public int numTriangles; 
	public boolean triangulated = false; 
	public int[][] gComponentsSingle = new int[6][6];
	public int[] cOldX = new int[3];
	public int[] cOldY = new int[3]; 
	public int[] cNewX = new int[3]; 
	public int[] cNewY = new int[3]; 
	public int numConstraints = 0; 

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
		assert(cOldX.length == cOldY.length);
		assert(cOldX.length == cNewX.length);
		assert(cNewX.length == cNewY.length);
		while(num > cOldX.length) {
			// Double array capacity
			int[] copyOldX = new int[cOldX.length*2]; 
			int[] copyOldY = new int[cOldX.length*2];
			int[] copyNewX = new int[cOldX.length*2];
			int[] copyNewY = new int[cOldX.length*2];
			// Copy data across 
			System.arraycopy(cOldX, 0, copyOldX, 0, tempX.length);
			System.arraycopy(cOldY, 0, copyOldY, 0, tempX.length);
			System.arraycopy(cNewX, 0, copyNewX, 0, tempX.length);
			System.arraycopy(cNewY, 0, copyNewY, 0, tempY.length);
			
			// Reset references 
			cOldX = copyOldX;
			cOldY = copyOldY;
			cNewX = copyNewX;
			cNewY = copyNewY;
		}
		assert(cOldX.length == cOldY.length);
		assert(cOldX.length == cNewX.length);
		assert(cNewX.length == cNewY.length);
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

	public int[] addConstraint(int x, int y) {
		int pointX, pointY, distToPoint; 

		for(int i=0; i<numConstraints; i++) {
			distToPoint = (cOldX[i] - x)*(cOldX[i] - x) + (cOldY[i] - y)*(cOldY[i] - y); 
			if(distToPoint < 25) {
				return new int[]{0, 0, 0};
			}
		}
		for(int i=0; i<triangles.length; i++) {
			for(int j=0; j<triangles[i].npoints; j++){
				pointX = triangles[i].xpoints[j];
				pointY = triangles[i].ypoints[j];
				distToPoint = (pointX - x)*(pointX - x) + (pointY - y)*(pointY - y); 
				if(distToPoint < 25) {
					numConstraints++; 
					resizeConstraints(numConstraints);
					cOldX[numConstraints-1] = pointX;
					cOldY[numConstraints-1] = pointY;
					System.out.println("New Constraint: (" + pointX + ", " + pointY + ")");

					return new int[]{1, pointX, pointY};
				}
			}
		}
		return new int[]{0, 0, 0};
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
		// Iterate through triangles 
		// Compute v vector (Igarashi et al. Step One) 
		ArrayList<Integer> vListX = new ArrayList<Integer>(); 
		ArrayList<Integer> vListY = new ArrayList<Integer>();
		int x, y;
		boolean isConstraint, inV; 
		for(int i=0; i<numTriangles; i++) {
			for(int j=0; j<triangles[i].npoints; j++) {
				x = triangles[i].xpoints[j];
				y = triangles[i].ypoints[j];
				isConstraint = false;
				for(int k=0; k<numConstraints; k++) { 
					if(x == cOldX[k] && y == cOldY[k]) {
						isConstraint = true; 
						break;
					}
				}
				inV = false; 
				for(int k=0; k<vListX.size(); k++) {
					if(x == vListX.get(k) && y == vListY.get(k)) {
						inV = true; 
						break; 
					}
				}
				if(!inV) {
					if(isConstraint) {
						vListX.add(x);
						vListY.add(y);
					} else {
						vListX.add(0, x);
						vListY.add(0, y);
					}
				}
			}
		}
		assert(vListX.size() == vListY.size());
		int[] v = new int[vListX.size() * 2]; 
		for(int i=0; i<vListX.size(); i++){
			v[2*i] = vListX.get(i); 
			v[(2*i) + 1] = vListY.get(i);
		}

		// Build G Matrix 
		int[][] G = new int[v.length][v.length];
		int[] vIdx = new int[3]; 
		for(int i=0; i<numTriangles; i++) { 
			// Get triangle vertex indices in v 
			for(int j=0; j<triangles[i].npoints; j++) {
				for(int k=0; k<v.length; k+=2) {
					if(triangles[i].xpoints[j] == v[k] && triangles[i].ypoints[j] == v[k+1]) {
						vIdx[j] = k;
					}
				}
			}
			// Calculate associated G components 

		}
	
	}


}