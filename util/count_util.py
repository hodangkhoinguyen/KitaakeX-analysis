def isZeroOneAccepted(parsingLine: list[str], parentColumns: list[int] = []) -> bool:
    """
    Determine if the data line satisfies:
    + Parent columns only contain 0/0 or ./
    + Children columns has exactly one occurrence of 0/1 or 1/1 

    Args:
        parsingLine (_type_): _description_
        parentColumns (list, optional): _description_. Defaults to [].
    """
    # Check if parent columns only contain 0/0 or ./
    for parent in parentColumns:
        if not ((parsingLine[parent].startswith('0/0')) or (parsingLine[parent].startswith("./"))):
            return False

    # Check if the other columns only contains exactly 1 of 0/1 or 1/1
    count = 0
    for i in range(0, len(parsingLine)):
        if (parsingLine[i][0:3] == "0/1"):
            count += 1
        elif (parsingLine[i][0:3] == "1/1"):
            count += 1

        if (count > 1):
            return False

    if (count != 1):
        return False
    return True

def isNotContainZeroTwo(parsingLine: list[str]):
    """
    Determine whether parsing line doesn't contains any 0/2, 1/2, or 2/2
    """
    for i in parsingLine:
        if i.startswith("0/2") or (i.startswith("1/2")) or (i.startswith("2/2")):
            return False

    return True

def findIndexHit(parsingLine):
    """ 
    Find index hit by 0/1 or 1/1
    Warning: should run isZeroOneAccepted first to check if the row is legal
    Output: what type of hit (0/1 or 1/1) and the occuring index
    """
    for i in range(len(parsingLine)):
        if (parsingLine[i].startswith("0/1")):
            return "0/1", i
        elif (parsingLine[i].startswith("1/1")):
            return "1/1.", i

    return None
