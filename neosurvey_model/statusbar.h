static int statusbar = 0;

char rotate_statusbar()
{
	statusbar++;
	if ( statusbar == 4 )
		statusbar = 0;
	
	switch (statusbar)
	{
		case 0:
			return '-';
			break;
		case 1:
			return '\\';
			break;
		case 2:
			return '|';
			break;
		case 3:
			return '/';
			break;
	}
	
	return '-';
}
