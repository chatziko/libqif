/*! \class Guessing
 *  \brief The guessing model of entropy.
 *
 *  For most information about this theory see
 */
class Guessing : public LeakageMeasure<double> {
	public:
		using LeakageMeasure::LeakageMeasure;

		double entropy(const prob& pi);

		double cond_entropy(const prob& pi);

		virtual const char* class_name() {
			return "Guessing";
		}
};
