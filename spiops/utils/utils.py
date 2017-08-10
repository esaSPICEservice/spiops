def valid_url(html_file_name):
    """
    This function returns a valid URL for an HTML given a filename. 
    The filename is checked in such way that URL non valid characters
    are replaced by other characters. 
    
    This was used due to the fact that we were using the following string:
    
       '67P/CG' -> '67P-CG'
    
    as part of an URL for 2D plotting and it would not work
    
    :param html_file_name: Input filename 
    :type html_file_name: str
    :return: Corrected Input filename without URL non valid characters
    :rtype: str
    """""

    for element in ['$', '_', '.', '+', '!', '*', '(', ')', '/', '\\']:
        
        if element in ['/', '\\', '!', '*', '$']:
            replacement = '-'
        else:
            replacement = '_'
        
        html_file_name = html_file_name.replace(element,replacement)

    return html_file_name


def convert_ESOCorbit2list():

    return


def convert_OEM2list():

    return