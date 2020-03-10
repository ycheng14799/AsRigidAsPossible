import java.awt.*;
import java.awt.event.*;
import java.applet.Applet;
import java.io.*;

/* Class AsRigidAsPossible */ 
public class AsRigidAsPossible extends Applet {
	static Frame myFrame = null; 
	Button bTriangulate; 
	Button bClear; 
	Button bQuit; 
	Panel mainPanel; 
	private MyGraphics myG; 
	public int borderSize; 
	public Color paintColor, bkColor; 

	public AsRigidAsPossible() { // Constructor
		borderSize = 10000;
        paintColor = Color.red;
        bkColor    = Color.lightGray;
		bTriangulate = new Button("Triangulate");
		bClear = new Button("Clear");
	}

	public void clearMe() {
		Graphics g = getGraphics(); 
		Dimension dimension = getSize();
		g.setColor(getBackground());
		g.fillRect(0, 0, dimension.width, dimension.height);
	}

	public void init() { // Initialize the applet 
		setBackground(Color.gray);
		setLayout(new BorderLayout());
		
		// Main panel 
		mainPanel = new Panel(); 
			mainPanel.setBackground(Color.lightGray);
		
		// UI
		Panel panel2 = new Panel(); 
			panel2.setBackground(Color.lightGray);
		bTriangulate = new Button("Triangulate");
			panel2.add(bTriangulate);
		bClear = new Button("Clear"); 
			panel2.add(bClear);
		if (myFrame != null) {
			bQuit = new Button("Quit");
			panel2.add(bQuit);
		}

		add("North", panel2);
		add("South", mainPanel);

		myG = new MyGraphics(paintColor, Color.yellow, Color.gray);
		
		// Triangulation function
		bTriangulate.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				Graphics g = getGraphics();
				myG.triangulate(g);
				g.setPaintMode();
			}
		});

		// Clear function 
		bClear.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				// Clear components 
				myG.clear(mainPanel, getBackground()); 
				// Clear screen 
				clearMe(); 
				Graphics g = getGraphics(); 
				g.setPaintMode();
			}
		});

		// Quit function 
		if (myFrame != null) {
		  bQuit.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				   System.exit(0);
			}
		  });
		}

		// Handling mouse input 
		this.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				Graphics g = getGraphics(); 
				g.setColor(paintColor);
				int x = e.getX();
				int y = e.getY();
				// Debug: Print mouse input coordinates 
				// System.out.println("(" + x + ", " + y + ")");
				myG.addVertex(g, x, y);
				g.setPaintMode(); 
			}
		});
	}
	
	public static void main(String[] args) { // Main
		myFrame = new Frame("As-Rigid-As-Possible Shape Manipulation Application"); // Create frame with title
		AsRigidAsPossible myApp = new AsRigidAsPossible(); // Initialize application
		myApp.init(); 

		// Add applet to the frame
		myFrame.add(myApp, BorderLayout.CENTER); // Add applet to frame 
		myFrame.pack(); // Set window size 
		myFrame.setSize(600, 500);
		myFrame.setVisible(true); // Make frame visible 

	} 
}