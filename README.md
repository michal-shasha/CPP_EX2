                                       
//325763498
//michalshasha8@gmail.com

                                       **מטלה 2- אופרטורים**
                                                                     
                                                                      הסבר על האופרטורים שמימשתי:
**operator+(Graph &other)-**  פונקציה זו מחברת שני גרפים.
 קודם היא בודקת אם המטריצות הן באותו גודל. אם לא, נזרקת שגיאה. 
 אחרת, היא יוצרת מטריצה חדשה שבה כל איבר הוא סכום האיברים המקבילים בשתי המטריצות של הגרפים.
 לאחר מכן, היא מאפסת את האלכסון הראשי לאפס (כדי לוודא שזה גרף חוקי). 
**operator-(Graph &other)-**פונקציה זו מחסרת שני גרפים.
 כמו בחיבור, היא בודקת אם המטריצות באותו גודל, ואם לא, נזרקת שגיאה. 
 אחרת, היא יוצרת מטריצה חדשה שבה כל איבר הוא ההפרש בין האיברים המקבילים במטריצות של הגרפים.
 גם כאן, האלכסון הראשי מאופס. 
**operator+()-** חיבור אונרי, פונקציה זו מחזירה את הגרף כפי שהוא (לא משנה דבר).     
**operator-()-** חיסור אונרי, פונקציה זו מחזירה גרף שבו כל המשקלים של המטריצה שונו לערכם השלילי.  
**operator+=(int scalar)-** פונקציה זו מוסיפה מספר קבוע לכל איבר במטריצת הסמיכות של הגרף. האלכסון הראשי מאופס לאחר החיבור.
**operator-=(int scalar)-**פונקציה זו מחסרת מספר קבוע מכל איבר במטריצת הסמיכות של הגרף. האלכסון הראשי מאופס לאחר החיסור.
**operator*=(int scalar)-**פונקציה זו מכפילה כל איבר במטריצת הסמיכות של הגרף במספר קבוע.
**operator/=(int scalar)-**פונקציה זו מחלקת כל איבר במטריצת הסמיכות של הגרף במספר קבוע. אם המספר הוא אפס, נזרקת שגיאה.
**operator*(Graph &other)-**פונקציה זו מכפילה שני גרפים (מכפלת מטריצות סמיכות). היא בודקת אם מספר העמודות במטריצה הראשונה שווה למספר השורות במטריצה השנייה, אחרת נזרקת שגיאה. לאחר מכן, היא מבצעת את הכפלת המטריצות ומאפסת את האלכסון הראשי.
**operator==(Graph &other)-**פונקציה זו בודקת אם שני גרפים שווים. היא בודקת את מספר הקודקודים ואת ערכי המטריצות. אם יש הבדל, היא מחזירה false.
**operator!=(Graph &other)-**פונקציה זו מחזירה את ההפך מהתוצאה של השוואת השוויון (==).
**isGraphContained(Graph &smallerGraph, Graph &largerGraph)-**פונקציה זו בודקת אם גרף קטן מוכל בגרף גדול על ידי בדיקת אם כל המשקלים במטריצת הסמיכות של הגרף הקטן תואמים לאלו שבגרף הגדול.
**operator>(Graph &other)-** פונקציה זו בודקת אם גרף אחד גדול מגרף אחר לפי סדר הכלה. אם אין הכלה, היא משווה את מספר הצלעות ואם מספר הצלעות שווה היא משווה את גודל המטריצה.
**operator<=(Graph &other)-**פונקציה זו בודקת אם גרף אחד קטן או שווה לאחר, לפי סדר הכלה ואם הם שווים אז לפי מספר הצלעות ואם מספר הצלעות שווה אז משווים לפי גודל המטריצה.
**operator<(Graph &other)-**פונקציה זו בודקת אם גרף אחד קטן מגרף אחר לפי סדר הכלה. אם אין הכלה, היא משווה את מספר הצלעות ואם יש צורך את גודל המטריצה.
**operator>=(Graph &other)-**פונקציה זו בודקת אם גרף אחד גדול או שווה לאחר, לפי סדר הכלה מספר הצלעות והשואת גודל המטריצות.
**operator++()-**פונקציה זו מוסיפה 1 לכל איבר במטריצת הסמיכות של הגרף. מכיוון שזהו אופרטור פרה-פיקס, ההוספה מתבצעת לפני שהערך מוחזר. הפונקציה מחזירה את האובייקט עצמו (*this) לאחר ביצוע ההוספה.
**operator--()-**פונקציה זו מפחיתה 1לכל איבר במטריצת הסמיכות של הגרף. כמו באופרטור הפרה-פיקס להוספה, ההפחתה מתבצעת לפני שהערך מוחזר. הפונקציה מחזירה את האובייקט עצמו (*this) לאחר ביצוע ההפחתה.
**operator++(int)**פונקציה זו מוסיפה 1 לכל איבר במטריצת הסמיכות של הגרף, אבל ההוספה מתבצעת אחרי שהערך המקורי מוחזר. הפונקציה מחזירה עותק של הגרף כפי שהיה לפני ביצוע ההוספה. לכן, הפונקציה מקבלת פרמטר int דמה המבדיל בינה לבין אופרטור הפרה-פיקס.
**operator--(int)-**פונקציה זו מפחיתה 1 לכל איבר במטריצת הסמיכות של הגרף, אבל ההפחתה מתבצעת אחרי שהערך המקורי מוחזר. הפונקציה מחזירה עותק של הגרף כפי שהיה לפני ביצוע ההפחתה. כמו באופרטור הפוסט-פיקס להוספה, הפונקציה מקבלת פרמטר int דמה המבדיל בינה לבין אופרטור הפרה-פיקס.
**:ostream &operator<<(std::ostream &out, Graph &graph)-**פונקציה זו מדפיסה את מטריצת הסמיכות של הגרף בצורה פורמטית.