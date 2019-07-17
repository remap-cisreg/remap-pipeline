while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        echo "$package - Take ass and kick name !!"
                        echo " "
                        echo "$package [options] application [arguments]"
                        echo " "
                        echo "options:"
                        echo "-h, --help                show brief help"
                        echo "-wd, --working-directory=WORK_DIR       specify a working directory where ReMap is"
                        echo "-fd, --file-modification=FILE_DIFF      specify a tab difference file : old names   new_name"
                        echo "-fm, --file-metadata=FILE_META          specify tab separated metadata file"
                        echo "-t, --tab-dir=DIR_TAB                   specify tab directory name"
                        exit 0
                        ;;
                -wd)
                        shift
                        if test $# -gt 0; then
                                export WORKING_DIR=$1
                        else
                                echo "no working directory specified"
                                exit 1
                        fi
                        shift
                        ;;
                --working-directory*)
                        export WORKING_DIR=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -fd)
                        shift
                        if test $# -gt 0; then
                                export FILE_DIFF=$1
                        else
                                echo "no modification file specified"
                                exit 1
                        fi
                        shift
                        ;;
                --file-modification*)
                        export FILE_DIFF=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -fm)
                        shift
                        if test $# -gt 0; then
                                export FILE_META=$1
                        else
                                echo "no modification file specified"
                                exit 1
                        fi
                        shift
                        ;;
                --file-metadata*)
                        export FILE_META=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -t)
                        shift
                        if test $# -gt 0; then
                                export DIR_TAB=$1
                        else
                                echo "no tab dir specified"
                                exit 1
                        fi
                        shift
                        ;;
                --tab-dir*)
                        export DIR_TAB=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                *)
                        break
                        ;;
        esac
done



# Checking if all argument are here
if [ -z "$WORKING_DIR" ]
then
      echo "missing working directory"
      exit 1
fi

if [ -z "$FILE_DIFF" ]
then
      echo "missing modification file"
      exit 1
fi

# if [ -z "$FILE_META" ]
# then
#       echo "missing metadata file"
#       exit 1
# fi


# echo $FILE_DIFF
# echo $WORKING_DIR
# echo $FILE_META

# Uniformize path to work dir with / at the end
if [[ $WORKING_DIR = */ ]]
then
	PATH_WORKING_DIR="$WORKING_DIR"
else
  PATH_WORKING_DIR = "$WORKING_DIR/"
fi

# Create absolute path to modification file
if [[ $FILE_DIFF = /* ]]
then
	PATH_FILE_DIFF = "$FILE_DIFF"
else
	PATH_FILE_DIFF="$PATH_WORKING_DIR$FILE_DIFF"
fi

#----------------------------------------------
# Dealing metadata file
#----------------------------------------------


# Create absolute path to metadata file
if [[ $FILE_META = /* ]]
then
	PATH_FILE_META = "$FILE_META"
else
	PATH_FILE_META="$PATH_WORKING_DIR$FILE_META"
fi

echo $PATH_WORKING_DIR
echo $PATH_FILE_DIFF
echo $PATH_FILE_META

# Check if files exists
if [ ! -f $PATH_FILE_DIFF ]
then
    echo "$PATH_FILE_DIFF file not found!"
    exit
fi

# if [ ! -f $PATH_FILE_META ]
# then
#   echo "$PATH_FILE_META file not found!"
#   exit
# fi


#----------------------------------------------
# Dealing with tab
#----------------------------------------------

if [[ $DIR_TAB = */ ]]
then
	TEMP_PATH_DIR_TAB="$DIR_TAB"
else
  TEMP_PATH_DIR_TAB="$DIR_TAB/"
fi


if [[ $TEMP_PATH_DIR_TAB == *"1.metadata"* ]]
then
  PATH_DIR_TAB="$PATH_WORKING_DIR$TEMP_PATH_DIR_TAB"
else
  PATH_DIR_TAB=$PATH_WORKING_DIR"1.metadata/"$TEMP_PATH_DIR_TAB
fi

# echo $PATH_DIR_TAB
find $PATH_DIR_TAB -name "*_summary.tab"
