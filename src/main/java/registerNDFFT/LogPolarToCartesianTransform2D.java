package registerNDFFT;

import net.imglib2.RealLocalizable;
import net.imglib2.RealPositionable;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.InvertibleRealTransform;
import net.imglib2.realtransform.PolarToCartesianTransform2D;

public class LogPolarToCartesianTransform2D implements InvertibleRealTransform
{
	private static double x( final double r, final double t )
	{
		return Math.exp(r) * Math.cos( t );
	}

	private static double y( final double r, final double t )
	{
		return Math.exp(r) * Math.sin( t );
	}

	private static double r( final double x, final double y )
	{
		return Math.log(Math.sqrt( x * x + y * y ));
	}

	private static double t( final double x, final double y )
	{
		return Math.atan2( y, x );
	}

	private final InverseRealTransform inverse;

	public LogPolarToCartesianTransform2D()
	{
		inverse = new InverseRealTransform( this );
	}

	@Override
	public int numSourceDimensions()
	{
		return 2;
	}

	@Override
	public int numTargetDimensions()
	{
		return 2;
	}

	@Override
	public void apply( final double[] source, final double[] target )
	{
		final double r = source[ 0 ];
		final double t = source[ 1 ];
		target[ 0 ] = x( r, t );
		target[ 1 ] = y( r, t );
	}

	@Override
	public void apply( final float[] source, final float[] target )
	{
		final double r = source[ 0 ];
		final double t = source[ 1 ];
		target[ 0 ] = ( float ) x( r, t );
		target[ 1 ] = ( float ) y( r, t );
	}

	@Override
	public void apply( final RealLocalizable source, final RealPositionable target )
	{
		final double r = source.getDoublePosition( 0 );
		final double t = source.getDoublePosition( 1 );
		target.setPosition( x( r, t ), 0 );
		target.setPosition( y( r, t ), 1 );
	}

	@Override
	public void applyInverse( final double[] source, final double[] target )
	{
		final double x = target[ 0 ];
		final double y = target[ 1 ];
		source[ 0 ] = r( x, y );
		source[ 1 ] = t( x, y );
	}

	@Override
	public void applyInverse( final float[] source, final float[] target )
	{
		final double x = target[ 0 ];
		final double y = target[ 1 ];
		source[ 0 ] = ( float ) r( x, y );
		source[ 1 ] = ( float ) t( x, y );
	}

	@Override
	public void applyInverse( final RealPositionable source, final RealLocalizable target )
	{
		final double x = target.getDoublePosition( 0 );
		final double y = target.getDoublePosition( 1 );
		source.setPosition( r( x, y ), 0 );
		source.setPosition( t( x, y ), 1 );
	}

	@Override
	public InvertibleRealTransform inverse()
	{
		return inverse;
	}

	@Override
	public LogPolarToCartesianTransform2D copy()
	{
		return this;
	}
}


