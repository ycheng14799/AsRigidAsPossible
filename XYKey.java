public class XYKey {
	private final int x; 
	private final int y; 

	public XYKey(int x, int y) {
		this.x = x;
		this.y = y; 
	} 

	@Override 
	public boolean equals(Object o) {
		if (this == o) return true;
        if (!(o instanceof XYKey)) return false;
        XYKey key = (XYKey) o;
        return x == key.x && y == key.y;
    }

	@Override
    public int hashCode() {
        int result = x;
        result = 31 * result + y;
        return result;
    }
}