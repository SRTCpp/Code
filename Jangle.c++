#include "Jangle.h"


latangle::latangle(angle a, latconvention l)
	: _lat(a), _convention(l)
{}	

latangle::latangle(const latangle& a)
		: _lat(a.north_positive()), _convention(NORTH_POSITIVE) 
{} 

angle latangle::north_positive() const
{
	angle answer(_lat);
	if (_convention == NORTH_NEGATIVE) answer*=-1;
	return answer;
}

angle latangle::north_negative() const
{
	angle answer(_lat);
	if (_convention == NORTH_POSITIVE) answer*=-1;
	return answer;
}

latangle latangle::as(latconvention l) const
{
	latangle answer(*this);
	if (l==NORTH_POSITIVE) answer = latangle(north_positive(),l);
	else if (l==NORTH_NEGATIVE) answer = latangle(north_negative(),l);
	return answer;
}

latangle& latangle::to(latconvention l)
{
	return (*this)=as(l);
}


lonangle::lonangle(angle a, lonconvention l)
		: _lon(a), _convention(l)
{}


lonangle::lonangle(const lonangle &a)
		: _lon(a._lon), _convention(a._convention) 
{}


angle lonangle::center_negative180() const
{
	//std::cout <<"Calling a lonangle to output in convention: CENTER_NEGATIVE180\n";
	
	angle answer(_lon);
	if (_convention == CENTER_0) {
		if (answer.degrees() >0. ) answer=angle(answer.degrees()-360.,angle::DEG);
	}

	else if (_convention ==CENTER_POSITIVE180) answer= -answer;
	return answer;
}

angle lonangle::center_positive180() const
{
	//std::cout <<"Calling a lonangle to output in convention: CENTER_POSITIVE180\n";
	
	angle answer(_lon);
	if (_convention == CENTER_0) {
		if (answer.degrees() > 0. ) answer = angle(360.-answer.degrees(),angle::DEG);
		else if (answer.degrees() < 0. ) answer = -answer;
	}
	else if (_convention == CENTER_NEGATIVE180 && answer.degrees() !=0.) answer= -answer;	
	return answer;
}


angle lonangle::center_0() const
{
	
	angle answer(_lon);
	if (_convention == CENTER_NEGATIVE180) {
		//std::cout <<answer.degrees();
		if (answer.degrees() < -180.) answer=angle(answer.degrees()+360.,angle::DEG);
		//std::cout <<" -> "<<answer.degrees()<<"\n";
	}
	else if(_convention == CENTER_POSITIVE180){
		if (answer.degrees() > 180.) answer=angle(360.-answer.degrees(),angle::DEG);
	}
	return answer;
}


lonangle lonangle::as(lonconvention l) const
{
	lonangle answer(*this);
	if (l==CENTER_0) answer= lonangle(center_negative180(),l);
	else if (l==CENTER_NEGATIVE180) answer=lonangle(center_0(),l);
	else if (l == CENTER_POSITIVE180) answer = lonangle(center_positive180(),l);
	return answer;
}

lonangle& lonangle::to(lonconvention l) 
{
	return (*this)=as(l);
}
